// Authors: Gilles Civario and Michael Lysaght
// Copyright 2014 Irish Centre for High-End Computing under MIT license
// See license.txt for the full license of this software
//
// Second step of improvement: the code now will only send to the default device
// a fair share of work, keeping some for processing on the host. The load balancing
// is dynamic and based on the time taken to process the share of work on each sides
// during the previous iteration of the loop.

#include <iostream>
#include <cstdlib>
#include <cmath>

#include <omp.h>

using namespace std;

#ifdef DOUBLE
    typedef double real;
    #define Sqrt( x ) sqrt( x )
#else
    typedef float real;
    #define Sqrt( x ) sqrtf( x )
#endif

#pragma omp declare target
const real G = 6.67384e-11;

real *x, *y, *z, *vx, *vy, *vz, *m;

void randomFill( real *&x, real inf, real sup, size_t n ) {
    x = new real[n];
    const real mult = (sup - inf) / real( RAND_MAX );
    for ( size_t i = 0; i != n; ++i ) {
        x[i] = rand() * mult + inf;
    }
}

void bodies( size_t n ) {
    randomFill( x, -100, 100, n );
    randomFill( y, -100, 100, n );
    randomFill( z, -100, 100, n );
    randomFill( vx, -10, 10, n );
    randomFill( vy, -10, 10, n );
    randomFill( vz, -10, 10, n );
    randomFill( m, 1, 10000, n );    
}

void cleanBodies() {
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] vx;
    delete[] vy;
    delete[] vz;
    delete[] m;
}

// Updated interface for the Newton function: it now takes two extra parameters
// n0 and n1 corresponding respectively to the first and last indexes to compute 
void Newton( size_t n0, size_t n1, size_t n, real dt ) {
    const real dtG = dt * G;
    // Now the work is done on the device if and only if the first index to
    // process is 0. The work is done on the host otherwise
    #pragma omp target if( n0 == 0 )
    #pragma omp parallel
    {
        #pragma omp for schedule( auto )
        for ( size_t i = n0; i < n1; ++i ) {
            real dvx = 0, dvy = 0, dvz = 0;
            #pragma omp simd
            for ( size_t j = 0; j < i; ++j ) {
                real dx = x[j] - x[i], dy = y[j] - y[i], dz = z[j] - z[i];
                real dist2 = dx*dx + dy*dy + dz*dz;
                real mOverDist3 = m[j] / (dist2 * Sqrt( dist2 ));
                dvx += mOverDist3 * dx;
                dvy += mOverDist3 * dy;
                dvz += mOverDist3 * dz;
            }
            #pragma omp simd
            for ( size_t j = i+1; j < n; ++j ) {
                real dx = x[j] - x[i], dy = y[j] - y[i], dz = z[j] - z[i];
                real dist2 = dx*dx + dy*dy + dz*dz;
                real mOverDist3 = m[j] / (dist2 * Sqrt( dist2 ));
                dvx += mOverDist3 * dx;
                dvy += mOverDist3 * dy;
                dvz += mOverDist3 * dz;
            }
            vx[i] += dvx * dtG;
            vy[i] += dvy * dtG;
            vz[i] += dvz * dtG;
        }
        #pragma omp for simd schedule( auto )
        for ( size_t i = n0; i < n1; ++i ) {
            x[i] += vx[i] * dt;
            y[i] += vy[i] * dt;
            z[i] += vz[i] * dt;
        }
    }
}
#pragma omp end declare target

int main( int argc, char *argv[] ) {
    const size_t n = 50000;

    bodies( n );

    // Force the possibility to have nested levels of parallelism
    // This is necessary since we have two nested OpenMP parallel sections: one
    // in the main function (a few lines down from here), and one inside the Newton
    // function. Without this call (or if the environment variable OMP_NESTED isn't
    // set properly at run-time), only the first parallel section would spawn new threads
    omp_set_nested( true );

    // Ratio of data to send to the device, initially set to 0.5
    double ratio = 0.5;

    // Table storing the compute time on the device and on the host
    double tth[2];

    cout << "Before to start:\n"
         << "  Position of first particle is (" << x[0] << ',' << y[0] << ',' << z[0] << ")\n"
         << "  Position of last particle is (" << x[n-1] << ',' << y[n-1] << ',' << z[n-1] << ")\n";

    double tm = omp_get_wtime();

    // Sending to the device all the necessary inputs. The retrieval of outputs will be handled
    // inside the loop itself
    #pragma omp target data map( to: x[0:n], y[0:n], z[0:n], vx[0:n], vy[0:n], vz[0:n], m[0:n] )
    {
        // First level of parallelism: spawning only 2 threads, the first to "manage" the device
        // and the second to "manage" the host
        #pragma omp parallel num_threads( 2 )
        {
            // Retrieving the thread index: tid 0, will offload to the device and tid 1 will keep
            // its share on the host 
            const int tid = omp_get_thread_num();
            for ( int it = 0; it < 100; ++it ) {
                // Threshold index between what goes on the device and what stays on the host
                size_t lim = n * ratio;
                // Number of elements to keep on the device
                size_t l = n - lim;

                double tt = omp_get_wtime();
                // Calling our new Newton function with the starting and ending indexes as first parameters
                Newton( lim*tid, lim + l*tid, n, 0.01 );
                tth[tid] = omp_get_wtime() - tt;

                // Need a barrier here to ensure that both sides have finish computing prior to start
                // synchronising the data between them
                #pragma omp barrier
                // Only one thread do the data exchanges (to avoid data collisions) while the other one waits
                // for completion at the end of the section
                #pragma omp single
                {
                    // Retrieving on the host the data computed on the device 
                    #pragma omp target update from( x[0:lim], y[0:lim], z[0:lim], vx[0:lim], vy[0:lim], vz[0:lim]  )
                    // Sending to the host the data computed on the host
                    // ATTENTION: notation for data transfers is of the form array[first_index:length_of_the_array]
                    #pragma omp target update to( x[lim:l], y[lim:l], z[lim:l], vx[lim:l], vy[lim:l], vz[lim:l] )
                    // Computing of the new ratio of data to offload, based on the measured compute times
                    ratio = ratio*tth[1] / ((1-ratio)*tth[0] + ratio*tth[1]);
                }
                // There is an implicit barrier here
            }
        }
    }
    tm = omp_get_wtime() - tm;

    cout << "At the end:\n"
         << "  Position of first particle is (" << x[0] << ',' << y[0] << ',' << z[0] << ")\n"
         << "  Position of last particle is (" << x[n-1] << ',' << y[n-1] << ',' << z[n-1] << ")\n";
    cout << "Time was " << tm << "s\n";

    cleanBodies();

    return 0;
}
