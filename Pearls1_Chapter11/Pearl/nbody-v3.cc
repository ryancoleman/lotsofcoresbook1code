// Authors: Gilles Civario and Michael Lysaght
// Copyright 2014 Irish Centre for High-End Computing under MIT license
// See license.txt for the full license of this software
//
// Final step of improvement: now, we will use as many devices as we can and
// offload their fair share of work to them. The host will dedicate one management
// thread per device, and use all the other available threads to compute its own
// share of the work. Data transfers are done as scarcely as possible (only what needs
// to be updated gets updated) but not bi-directionally conversely to version 2.1. The
// data retrievals from the devices are done sequentially (one device at a time), but
// sends to all the devices are done in parallel.

#include <iostream>
#include <cstdlib>
#include <cmath>

#include <omp.h>

using namespace std;

#ifdef DOUBLE
    typedef double real;
    #define Sqrt sqrt
#else
    typedef float real;
    #define Sqrt sqrtf
#endif

#pragma omp declare target
real *x, *y, *z, *vx, *vy, *vz, *m;
const real G = 6.67384e-11;

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

// Yet again a modified version of the Newton function interface. Now, the parameters
// correspond to:
// n: size of the entire arrays (as usual)
// s: starting index of the sub-arrays to compute (somewhat new)
// l: length of the sub-arrays to compute (new)
// dt: time increment (as usual)
// tid: index of the calling thread (new)
// dev: number of devices to offload to (new)
void Newton( size_t n, size_t s, size_t l, real dt, int tid, int dev ) {
    const real dtG = dt * G;
    // Only offload some work to the device numbered "tid" if my thread index is lesser than the number
    // of devices "dev". Thus, thread #0 manages device #0, etc, and thread #dev manages the host
    #pragma omp target device( tid ) if( tid < dev )
    #pragma omp parallel
    {
        #pragma omp for schedule( auto )
        for ( size_t i = s; i < s+l; ++i ) {
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
        for ( size_t i = s; i < s+l; ++i ) {
            x[i] += vx[i] * dt;
            y[i] += vy[i] * dt;
            z[i] += vz[i] * dt;
        }
    }
}
#pragma omp end declare target

// New function to compute index displacements into the arrays for each of the devices
// (host is seen as the last device) depending on the times measured after the
// previous iteration. For a device of index d, the work that will be managed spans
// from index displ[d] to index displ[d+1]
inline void computeDisplacements( size_t *displ, const double *tth, int dev ) {
    double sumLengthOverT = 0;
    size_t length[dev+1];
    for ( int i = 0; i < dev+1; ++i ) {
        length[i] = displ[i+1] - displ[i];
    }
    for ( int i = 0; i < dev+1; ++i ) {
        sumLengthOverT += length[i] / tth[i];
    }
    for ( int i = 0; i < dev; ++i ) {
        displ[i+1] = displ[i] + round( (displ[dev+1] * length[i]) / (tth[i]*sumLengthOverT) );
    }
}

int main( int argc, char *argv[] ) {
    // Retrieving the number of devices available on the machine
    int dev = omp_get_num_devices();

    const size_t n = 50000;

    // Computing the initial index displacements into the arrays for each device
    // For a start, each device and the host get an equal share of the work
    size_t displ[dev+2];
    displ[0] = 0;
    for ( int i = 1; i < dev+1; ++i ) {
        displ[i] = (i * n) / (dev + 1);
    }
    displ[dev+1] = n;

    bodies( n );

    omp_set_nested( true );

    cout << "Before to start:\n"
         << "  Position of first particle is (" << x[0] << ',' << y[0] << ',' << z[0] << ")\n"
         << "  Position of last particle is (" << x[n-1] << ',' << y[n-1] << ',' << z[n-1] << ")\n";

    double tth[dev+1];
    double tm = omp_get_wtime();
    
    // Now, devices are managed in parallel, which means that, unless we explicitly request otherwise
    // (through a critical section for example), all data transfers are done in parallel
    #pragma omp parallel num_threads( dev + 1 )
    {
        const int tid = omp_get_thread_num();
        #pragma omp target data device( tid ) if( tid < dev ) \
            map( to: x[0:n], y[0:n], z[0:n], vx[0:n], vy[0:n], vz[0:n], m[0:n] )
        {
            for ( int it = 0; it < 100; ++it ) {
                // Computing the starting indexes and lengths of work to carry out on the current "device"
                // (this can be a device or the host)
                size_t s = displ[tid], l =  displ[tid+1] - displ[tid];

                // Doing the work
                double tt = omp_get_wtime();
                Newton( n, s, l, 0.01, tid, dev );
                tth[tid] = omp_get_wtime() - tt;

                // Now we retrieve sequentially on the host the various segments of data computed on the devices
                // (the sequential aspect comes from the critical section, and is for avoiding data collisions)
                #pragma omp critical
                {
                    #pragma omp target update device( tid ) if( tid < dev ) \
                        from( x[s:l], y[s:l], z[s:l], vx[s:l], vy[s:l], vz[s:l] )
                }

                // A barrier for data synchronisation
                #pragma omp barrier

                // Now on the host we have all the arrays up-to-date
                // We can therefore send back to the devices what is missing for them
                // which corresponds to at most two segments per array
                // We send the segment which is before the data the current device is responsible of
                #pragma omp target update device( tid ) if( tid < dev ) \
                    to( x[0:s], y[0:s], z[0:s], vx[0:s], vy[0:s], vz[0:s] )

                // then we compute the starting index and length of the segments following
                // the data the current device owns
                size_t s1 = s+l, l1 = n-s-l;

                // And we send this last missing segment
                #pragma omp target update device( tid ) if( tid < dev ) \
                    to( x[s1:l1], y[s1:l1], z[s1:l1], vx[s1:l1], vy[s1:l1], vz[s1:l1] )

                // Finally we compute the new shares of work to offload
                #pragma omp single
                computeDisplacements( displ, tth, dev );
                // And we have a final implicit barrier for data synchronisation here
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
