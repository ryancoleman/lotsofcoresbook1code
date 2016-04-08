#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <nlopt.h>
#include "mpi.h"

int numTasks, mpiRank;

#include "myFunc.h"

void writeParam(char *filename, int nParam, double *x)
{
  FILE *fn=fopen(filename,"w");
  if(!fn) {
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(1);
  }

  uint32_t n=nParam; // ensure size is uint32_t
  fwrite(&n,sizeof(uint32_t), 1, fn);
  for(int i=0; i < nParam; i++) {
    float tmp=x[i];
    fwrite(&tmp,sizeof(float), 1, fn);
  }
  fclose(fn);
}

int masterOP=1;
void startClient(void * restrict my_func_data)
{
  int op;
  double xFromMPI[N_PARAM];
  double partialError,sum;
  
  for(;;) { // loop until the master says I am done - then exit
    MPI_Bcast(&op, 1, MPI_INT, 0, MPI_COMM_WORLD); // receive the op code
    if(op==0) { // we are done, normal exit
      break;
    }
    MPI_Bcast(xFromMPI, N_PARAM, MPI_DOUBLE, 0, MPI_COMM_WORLD); // receive the parameters
    partialError = objFunc(N_PARAM,  xFromMPI, NULL, my_func_data);
    MPI_Reduce(&partialError, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
}

double mpiObjFunc(unsigned n,  const double * restrict x, double * restrict grad,
		  void * restrict my_func_data)
{
  int op;
  double partialError, totalError=0.;

  MPI_Bcast(&masterOP, 1, MPI_INT, 0, MPI_COMM_WORLD); // Send the master op code
  MPI_Bcast((void*) x, N_PARAM, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Send the parameters
  partialError = objFunc(N_PARAM,  x, NULL, my_func_data);
  MPI_Reduce(&partialError, &totalError, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // get the totalError

  return(totalError);
}

int main(int argc, char* argv[])
{
  nlopt_opt opt;
  userData_t uData = {0};
  FILE *fout_mpi;

  if(argc < 4) {
    fprintf(stderr,"Use: datafile paramFile mpiTimingFile\n");
    return -1;
  }

  int ret = MPI_Init(&argc,&argv);
  if (ret != MPI_SUCCESS) {
    printf ("Error in MPI_Init()!\n");
    MPI_Abort(MPI_COMM_WORLD, ret);
  }
  
  MPI_Comm_size(MPI_COMM_WORLD,&numTasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

#if defined(NO_IO_EXAMPLES)
#pragma message "TEST case: Clients use random memory values"
    init_noIO(NO_IO_EXAMPLES,&uData); 
#else
#pragma message "TEST case: all clients perform file reads"
  { // for simplicity, append the mpiRank to the data filename
    char filename[256];
    sprintf(filename,"%s.%d",argv[1],mpiRank);
    //fprintf(stderr,"Loading %s into coprocessor %d\n", filename, MIC_DEV);
    init(filename,&uData); 
  }
#endif

  if(mpiRank > 0) {
    startClient(&uData);
  } else { // Master code

    { // open MPI results file for this size run
      char buf[256];
      sprintf(buf, "%s.%04d.txt",argv[3],numTasks);
      fout_mpi = fopen(buf,"w");
      if(!fout_mpi) {
	fprintf(stderr,"Cannot open file %s\n",buf);
	return -1;
      }
    }

    printf ("Number of tasks= %d My rank= %d, number clients %d\n", numTasks,mpiRank,numTasks);
#ifdef MPI_NUM_COPROC_PER_NODE
    printf ("Number of coprocessors per node= %d\n", MPI_NUM_COPROC_PER_NODE);
#endif
    printf("myFunc %s\n", desc);
    printf("nExamples %d\n", uData.nExamples);
    printf("Number Parameters %d\n", N_PARAM);
    
    opt = nlopt_create(NLOPT_LN_PRAXIS, N_PARAM); // algorithm and dimensionality
    // NOTE: alternative optimization methods ...
    //opt = nlopt_create(NLOPT_LN_NEWUOA, N_PARAM);
    //opt = nlopt_create(NLOPT_LN_COBYLA, N_PARAM);
    //opt = nlopt_create(NLOPT_LN_BOBYQA, N_PARAM);
    //opt = nlopt_create(NLOPT_LN_AUGLAG, N_PARAM);
    
    nlopt_set_min_objective(opt, mpiObjFunc, (void*) &uData);
#if defined(MAX_RUNTIME) 
  fprintf(stderr,"Warning: performing a quick %d second test!\n", MAX_RUNTIME);
  nlopt_set_maxtime(opt, MAX_RUNTIME); // Use for running quick tests
#else
  fprintf(stderr,"MAX runtime %d seconds!\n", 120*60);
  nlopt_set_maxtime(opt, (120. * 60.)); // maximum runtime in seconds
#endif

    double minf; /* the minimum objective value, upon return */
    
    __declspec(align(64)) double x[N_PARAM];
    for(int i=0; i < N_PARAM; i++) x[i] = 0.1*(random()/(double)RAND_MAX);
    
    double startTime=getTime();
    ret=nlopt_optimize(opt, x, &minf);
    printf("Optimization Time %g\n",getTime()-startTime);
    
    if (ret < 0) {
      printf("nlopt failed! ret %d\n", ret);
    } else {
      printf("found minimum %0.10g ret %d\n", minf,ret);
    }
    writeParam(argv[2],N_PARAM, x);

    nlopt_destroy(opt);
    
    masterOP = 0; // signal completion
    MPI_Bcast(&masterOP, 1, MPI_INT, 0, MPI_COMM_WORLD); // Send the master op code
    printf("----------- performance times for the Master ----------\n");
    fini(&uData);
  }

  int client_nExamples[numTasks];
  double client_timeObjFunc[numTasks];
  int client_countObjFunc[numTasks];
  double client_timeDataLoad[numTasks];
  double client_minTime[numTasks];
  double client_maxTime[numTasks];
  
  MPI_Gather(&uData.nExamples, 1, MPI_INT, client_nExamples, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&uData.countObjFunc, 1, MPI_INT, client_countObjFunc, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&uData.timeDataLoad, 1, MPI_DOUBLE, client_timeDataLoad, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&uData.timeObjFunc, 1, MPI_DOUBLE, client_timeObjFunc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&uData.minTime, 1, MPI_DOUBLE, client_minTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&uData.maxTime, 1, MPI_DOUBLE, client_maxTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(mpiRank==0) {
    printf("----------- performance times for the MPI run ----------\n");
    printf("function: %s\n",desc);
    
    uint64_t totalExamples=0;
    for(int i=0; i < numTasks; i++) totalExamples += client_nExamples[i];
    printf("totalExamples %g\n",(double) totalExamples);
    
    printf("AveObjTime %g, countObjFunc %d, totalObjTime %g\n",
	   uData.timeObjFunc/uData.countObjFunc, uData.countObjFunc, uData.timeObjFunc);
#ifdef FLOP_ESTIMATE
    printf("Estimated flops in myFunc %d, estimated average TFlop/s %g, nClients %d\n", FLOP_ESTIMATE,
	   (((double)totalExamples * (double)FLOP_ESTIMATE)/(uData.timeObjFunc/uData.countObjFunc)/1.e12),
	   numTasks);
    printf("Estimated maximum TFlop/s %g, minimum TFLop/s %g\n",
	   (((double)totalExamples*(double)FLOP_ESTIMATE)/(uData.minTime)/1.e12),
	   (((double)totalExamples*(double)FLOP_ESTIMATE)/(uData.maxTime)/1.e12) );
#endif
    
    fprintf(fout_mpi, "nExamples countObjFunc timeDataLoad  timeObjFunc minObjFcnTime maxObjFcnTime\n");
    for(int i=0; i < numTasks; i++) {
      fprintf(fout_mpi, "%g %g %g %g %g %g\n",
	      (double) client_nExamples[i], (double) client_countObjFunc[i], client_timeDataLoad[i],
	      client_timeObjFunc[i], client_minTime[i], client_maxTime[i]);
    }
  }

  MPI_Finalize();
  
  return 0;
}
