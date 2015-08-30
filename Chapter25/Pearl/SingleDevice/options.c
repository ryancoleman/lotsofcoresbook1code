/*
 * =====================================================================================
 *
 *       Filename:  options.cpp
 *
 *    Description:  Asian put/call options MC simulator
 *
 *        Version:  1.1
 *        Created:  
 *       Revision:  none
 *       Compiler:  icpc
 *
 *         Author:  Vadim Karpusenko, PhD (vk), vadikus@gmail.com
 *        Company:  Colfax International
 *
 * =====================================================================================
 */

#include <mpi.h>
#include <mkl.h>
#include <mkl_vsl.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <unistd.h>

// Constants for MPI communication
const int bossRank = 0;
const int msgReportLength = 8;
const int msgReportTag = 1;
const int msgSchedLength = 8;
const int msgSchedTag    = 2;
const int msgNameTag     = 3;
#define hostNameLen 128
typedef char HostNameType[128];

// Constants for
#define verboseCommunication 0
#define verboseTiming 1
#define verboseResults 1

#define ALIGNED __attribute__((aligned(64)));

int main(int argc, char **argv) {

  MPI_Status mpiStatus;
  int myRank, mpiWorldSize;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);

  if (myRank == bossRank) {

    HostNameType workerName[mpiWorldSize];
    for (int i = 1; i < mpiWorldSize; i++) {
      int len;
      MPI_Recv(workerName[i], hostNameLen, 
          MPI_CHAR, i, msgNameTag,
          MPI_COMM_WORLD, &mpiStatus);
    }

    const int nCalculations = 10;
    const int vMin = 0.05;
    const int vMax = 0.50;

    for (int iCalc = 0; iCalc < nCalculations; iCalc++) {

      const double t0 = omp_get_wtime();

      // Boss process
      const int nStrikes = 100;
      const float strikeMin = 10.0f;
      const float strikeMax = 20.0f;

      const int numPaths = 1<<20;
      const int numIntervals = 365;
      const float S = 15.3;
      const float r = 0.08;
      const float v = vMin + (vMax - vMin)*(float)(iCalc)/(float)(nCalculations - 1);
      const float T = 1.0;

      float payoff_geom_put   [nStrikes];
      float payoff_geom_call  [nStrikes];
      float payoff_arithm_put [nStrikes];
      float payoff_arithm_call[nStrikes];

      float workerTime[mpiWorldSize];
      float workerPath[mpiWorldSize];
      int   workerTask[mpiWorldSize];

      workerTime[:] = 0.0f;
      workerPath[:] = 0.0f;
      workerTask[:] = 0;

      // Loop over tasks, handing them out to 
      // available workers, one task at a time
      int nReceivedResults = 0;
      int iStrike = 0;
      while (nReceivedResults < nStrikes) {

        // Wait for any worker to report for work
        float msgBuf[msgReportLength];
        if (verboseCommunication) printf("Boss waiting for report...\n");
        MPI_Recv(&msgBuf, msgReportLength,
            MPI_INT, MPI_ANY_SOURCE, msgReportTag,
            MPI_COMM_WORLD, &mpiStatus);
        const int reportingWorker = mpiStatus.MPI_SOURCE;
        const int haveResults = floorf(msgBuf[0]);

        if (verboseCommunication)
          printf("Boss received report from worker %d %s\n",
              reportingWorker, (haveResults ? "(with results)" : ""));

        if (haveResults) {
          // Parse and record the results from reporing worker
          nReceivedResults++;

          const int iStrikeReported = floorf(msgBuf[1]); // The index number of the strike price

          workerTime[reportingWorker] += msgBuf[2]; // The amount of time the worker was computing
          workerPath[reportingWorker] += floorf(msgBuf[3]); // The number of paths processed by the worker
          workerTask[reportingWorker] += 1; // The number of tasks processed by the worker

          payoff_geom_put   [iStrikeReported] = msgBuf[4]; // Computed put  option payoff, geometric mean 
          payoff_geom_call  [iStrikeReported] = msgBuf[5]; // Computed call option payoff, geometric mean 
          payoff_arithm_put [iStrikeReported] = msgBuf[6]; // Computed put  option payoff, arithmetic mean
          payoff_arithm_call[iStrikeReported] = msgBuf[7]; // Computed call option payoff, arithmetic mean

        }

        if (iStrike < nStrikes) {
          // Assign the next task to the worker

          // Choose the next strike price
          const float K = strikeMin + (strikeMax - strikeMin)*
            (float)iStrike/(float)(nStrikes - 1);

          float msgBuf[msgSchedLength];
          msgBuf[0] = numPaths + 0.01f;     // Number of paths to simulate (0 to quit)
          msgBuf[1] = numIntervals + 0.01f; // Number of intervals for averaging
          msgBuf[2] = iStrike + 0.01f;      // Index number of the strike price
          msgBuf[3] = K;                    // Value of the strike price
          msgBuf[4] = S;		    // Starting asset value	 
          msgBuf[5] = r;		    // Risk-free rate		 
          msgBuf[6] = v;		    // Option volatility	 
          msgBuf[7] = T;		    // Option expiration time   

          // Sending floating-point parameters
          MPI_Send((void*)&msgBuf, msgSchedLength, MPI_FLOAT, reportingWorker, msgSchedTag, MPI_COMM_WORLD);
          if (verboseCommunication) printf("Boss sent command to worker %d (%d)...\n", reportingWorker, (int)floorf(msgBuf[0]));

          iStrike++;

        }
      }

      // Print calculation results
      if (verboseResults) {
        printf("Results:\n");
        printf("Number of paths: %d\n", numPaths);
        printf("Number of intervals: %d\n", numIntervals);
        printf("S=%.2f, r=%.3f, v=%.2f, T=%.2f\n", S, r, v, T);
        printf("%1s%8s%8s%8s%8s%8s%8s\n",
            "#", "i", "K", "G/Put", "G/Call", "A/Put", "A/Call");
        for (int iStrike = 0; iStrike < nStrikes; iStrike++) {

          const float K = strikeMin + (strikeMax - strikeMin)*
            (float)iStrike/(float)(nStrikes - 1);

          printf("%1s%8d%8.2f%8.2f%8.2f%8.2f%8.2f\n",
              " ", iStrike, K, 
              payoff_geom_put[iStrike], payoff_geom_call[iStrike],
              payoff_arithm_put[iStrike], payoff_arithm_call[iStrike]);
        }
      }

      const double t1 = omp_get_wtime();

      // Print timing results
      if (verboseTiming) {
        printf("%1s%15s%8s%12s%8s\n", "#", "Worker", "Share", "Performance", "Effic.");
        for (int i = 1; i < mpiWorldSize; i++) {
          printf("%16s%7.1f%%%12.2e%7.1f%%\n",
              workerName[i], 100.0f*(float)workerTask[i]/(float)nStrikes,
              (workerTime[i] > 0.0f ? workerPath[i]/workerTime[i] : 0.0f), 100.0f*workerTime[i]/(t1-t0));
        }
        printf("# Calculation %3d of %3d took %.3f seconds\n", iCalc+1, nCalculations, t1-t0);
        printf("# Net performance: %.2e paths/second\n\n", (float)numPaths*(float)nStrikes/(t1-t0));

        if (iCalc == nCalculations - 1) {
          printf("\n# Explanation of columns:\n\
# Worker      - the host on which the respective worker process was running\n\
# Share       - the fraction of work processed by that worker (%%)\n\
# Performance - the performance of the worker: number of option paths simulated per second\n\
# Effic.      - the fraction of time that the worker was loaded with calculations\n\n");
        }
      }


      for (int iWorker = 1; iWorker < mpiWorldSize; iWorker++) {
        // Telling workers to come back tomorrow (i.e., for the next iCalc)
        float msgBuf[msgSchedLength];
        msgBuf[0:msgSchedLength] = 0.0f;
        msgBuf[0] = 0.01f; // 0 paths is a signal to report to work again
        MPI_Send((void*)&msgBuf, msgSchedLength, MPI_INT, iWorker, msgSchedTag, MPI_COMM_WORLD);
        if (verboseCommunication) printf("Boss sent command to worker %d (%d)...\n", iWorker, (int)floorf(msgBuf[0]));
      }

    }

    for (int iWorker = 1; iWorker < mpiWorldSize; iWorker++) {
      // Terminating workers
      float msgBuf[msgSchedLength];
      msgBuf[0:msgSchedLength] = 0.0f;
      msgBuf[0] = -10.0f; // negative number of paths is a signal to quit
      MPI_Send((void*)&msgBuf, msgSchedLength, MPI_INT, iWorker, msgSchedTag, MPI_COMM_WORLD);
      if (verboseCommunication) printf("Boss sent command to worker %d (%d)...\n", iWorker, (int)floorf(msgBuf[0]));
    }
  
  // End Boss rank here
  
  } else {

    // Worker process code starts here

    char myName[hostNameLen];
    int len;
    MPI_Get_processor_name( &myName[0], &len);
    MPI_Send( (void*)&myName, hostNameLen, MPI_CHAR, bossRank,
        msgNameTag, MPI_COMM_WORLD);

    // Variables shared between threads:

    // ... parameters of calculation (input):
    int iStrike;
    int numPaths=1;         // Number of simulated asset paths
    int num_intervals;     // Number of intervals for the asset path to be sampled
    float S;               // Underlying asset price
    float K;               // Strike price
    float r;               // Risk-free rate
    float v;               // Volatility of the underlying asset
    float T;               // One year until expiration

    // ... payoffs of the option (output):
    float payoff_geom_put;    // geometric mean put options payoff
    float payoff_geom_call;   // geometric mean call options payoff
    float payoff_arithm_put;  // arithmetic mean put payoff
    float payoff_arithm_call; // arithmetic mean call payoff

    int haveResults = 0;
    int quit = 0;
    int set = 0;
    double execTime = 0.0;


    // Starting a thread-parallel region
    // which will exist until the boss terminates this worker
#pragma omp parallel
    {
      // For each thread, initialize a random number stream
      VSLStreamStatePtr stream; // stream for random numbers
      int seed = omp_get_thread_num();
      int errcode = vslNewStream( &stream, VSL_BRNG_MT2203, seed );

      // For each thread, initialize scratch variables
      // These quantities will evolve multiple random paths
      // in each thread using SIMD instructions (vector operations)
      const int vec_size = 1024;
      float spot_prices      [vec_size] ALIGNED;
      float rands            [vec_size] ALIGNED;
      float sumsm            [vec_size] ALIGNED;
      float sumsg            [vec_size] ALIGNED;
      float geom_mean        [vec_size] ALIGNED;
      float arithm_mean      [vec_size] ALIGNED;
      float geom_mean_put    [vec_size] ALIGNED;
      float arithm_mean_put  [vec_size] ALIGNED;
      float geom_mean_call   [vec_size] ALIGNED;
      float arithm_mean_call [vec_size] ALIGNED;

      while (numPaths >= 0) { // Repeat until the boss sends a command to quit

#pragma omp master
        {
          {
            // Report to boss: ask for work and
            // send back results of the previous calculation (if any),
            float msgBuf[msgReportLength];
            msgBuf[0:msgReportLength] = 0.0f;
            msgBuf[0] = haveResults + 0.1f;

            if (haveResults) {

              msgBuf[1] = iStrike + 0.1f;     // Index number of the strike price
              msgBuf[2] = execTime;           // The amount of time the worker was computing
              msgBuf[3] = (float)numPaths;    // The number of processed paths
              msgBuf[4] = payoff_geom_put;    // Computed put  option payoff, geometric mean
              msgBuf[5] = payoff_geom_call;   // Computed call option payoff, geometric mean
              msgBuf[6] = payoff_arithm_put;  // Computed put  option payoff, arithmetic mean
              msgBuf[7] = payoff_arithm_call; // Computed call option payoff, arithmetic mean

            }

            MPI_Send(&msgBuf, msgReportLength, MPI_FLOAT, bossRank, msgReportTag, MPI_COMM_WORLD);
          }

          {
            // Receive a scheduling command from boss
            float msgBuf[msgSchedLength];
            MPI_Recv(&msgBuf, msgSchedLength, MPI_FLOAT, bossRank, msgSchedTag, MPI_COMM_WORLD, &mpiStatus);

            numPaths      = floorf(msgBuf[0]); // Number of paths to simulate (0 to quit)
            num_intervals = floorf(msgBuf[1]); // Number of time intervals for averaging
            iStrike       = floorf(msgBuf[2]); // Index number of the strike price
            K             =       (msgBuf[3]); // Value of the strike price
            S             =       (msgBuf[4]); // Starting asset value
            r             =       (msgBuf[5]); // Risk-free rate
            v             =       (msgBuf[6]); // Option volatility
            T             =       (msgBuf[7]); // Option expiration time   

            if (numPaths > 0) execTime = omp_get_wtime();
            haveResults = 0;
          }
        }

        // A barrier is required so that non-master threads do not start
        // computing until the MPI command is received from the boss process
#pragma omp barrier 

        if (numPaths > 0) {

          // For each thread, initialize private variables for computation
          const int num_sims_vec = numPaths / vec_size; 
          const float dt = T / (float)num_intervals; // time interval between samples (in years)
          const float drift = dt*(r-0.5f*v*v);
          const float vol = v*sqrtf(dt);
          const float recipIntervals = 1.0f/(float)num_intervals;
          const float logS = logf(S);   

          // Team threads to process this loop in parallel 
          // The Asian option payoff calculation begins here
#pragma omp for schedule(guided)			\
          reduction(+: payoff_geom_put, payoff_geom_call,	\
              payoff_arithm_call, payoff_arithm_put)

          // Loop over the requested number of randomly simulated paths
          // (divided by the number of paths simultaneously computed in each thread
          // using vector instructions)
          for ( int i = 0; i < num_sims_vec; i++){

            // Each thread carries out calculations for multiple paths.
            // Initializing the state variables for these paths.
            for ( int k = 0; k < vec_size; k++) {
              spot_prices[k] = S; // initialize underlying assets
              sumsm[k] = S;       // initialize sum vector (arithmetic mean)
              sumsg[k] = logS;    // initialize sum vector (geometric mean)
            }

            // Loop over time intervals at which the Asian option price
            // is averaged
            for ( int j = 1; j < num_intervals; j++){
              vsRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, vec_size, rands, 0.0f, 1.0f);

              // Do the calculation for vec_size paths at in each thread
              // This loop is automatically vectorized by the compiler
              // This is loop is the performance-critical part of the calculation
              for ( int k = 0; k < vec_size; k++){

                spot_prices[k] *= expf(drift + vol*rands[k]); // Stochastic evolution of price
                sumsm[k] += spot_prices[k];         // arithmetic mean
                sumsg[k] += logf(spot_prices[k]);   // geometric mean

              }

            } // End of loop over time intervals

            // Computing the payoff for call and put options
            for ( int k = 0; k < vec_size; k++) {
              geom_mean_put[k]    = K - expf(sumsg[k] * recipIntervals); // put option
              geom_mean_call[k]   = - geom_mean_put[k];                  // call option
              arithm_mean_put[k]  = K - (sumsm[k] * recipIntervals);     // put option
              arithm_mean_call[k] = - arithm_mean_put[k];                // call option
              if (geom_mean_put[k]    < 0.0f) geom_mean_put[k]     = 0.0f;
              if (geom_mean_call[k]   < 0.0f) geom_mean_call[k]    = 0.0f;
              if (arithm_mean_call[k] < 0.0f) arithm_mean_call[k]  = 0.0f;
              if (arithm_mean_put[k]  < 0.0f) arithm_mean_put[k]   = 0.0f;
            }

            // Reduction of paths calculated in vector lanes
            // into scalar variables.
            // Simultaneously with reduction across vector lanes,
            // OpenMP reduces these scalars across threads
            for ( int k=0; k < vec_size; k++) {
              payoff_geom_put    += geom_mean_put[k];
              payoff_geom_call   += geom_mean_call[k];
              payoff_arithm_put  += arithm_mean_put[k];
              payoff_arithm_call += arithm_mean_call[k];
            }

          } // End of loop over the random paths

#pragma omp master
          {
            // Rescaling the quantities to convert sums into means
            payoff_geom_put    *= expf(-r*T)/((float)num_sims_vec*(float)vec_size);
            payoff_geom_call   *= expf(-r*T)/((float)num_sims_vec*(float)vec_size);
            payoff_arithm_put  *= expf(-r*T)/((float)num_sims_vec*(float)vec_size);
            payoff_arithm_call *= expf(-r*T)/((float)num_sims_vec*(float)vec_size);

            execTime = omp_get_wtime() - execTime;
            haveResults = 1; 
          }

          // End of Asian option payoff calculation
        }

#pragma omp barrier

      } // while(numPaths >= 0)

      vslDeleteStream( &stream );
    } // omp parallel

  } // worker process

  if (verboseCommunication) printf("Finishing. Rank : %d\n", myRank);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
