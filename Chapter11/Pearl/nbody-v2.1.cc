// Authors: Gilles Civario and Michael Lysaght
// Copyright 2014 Irish Centre for High-End Computing under MIT license
// See license.txt for the full license of this software
//
// Refinement of the second step of improvement: the data transfers between
// the host and the device are now done in parallel. This means that some data
// are sent to the device while some other are received on the host. This should
// exploit the bi-directional capability of the PCI express bus. However, since
// the computation times far dominate the data transfer times here, the expected
// performance improvement is minimal.

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

void Newton( size_t n0, size_t n1, size_t n, real dt ) {
    const real dtG = dt * G;
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

    omp_set_nested( true );

    double ratio = 0.5;
    double tth[2];

    cout << "Before to start:\n"
         << "  Position of first particle is (" << x[0] << ',' << y[0] << ',' << z[0] << ")\n"
         << "  Position of last particle is (" << x[n-1] << ',' << y[n-1] << ',' << z[n-1] << ")\n";

    double tm = omp_get_wtime();
    #pragma omp target data map( to: x[0:n], y[0:n], z[0:n], vx[0:n], vy[0:n], vz[0:n], m[0:n] )
    {
        #pragma omp parallel num_threads( 2 )
        {
            const int tid = omp_get_thread_num();
            for ( int it = 0; it < 100; ++it ) {
                size_t lim = n * ratio;
                size_t l = n - lim;

                double tt = omp_get_wtime();
                Newton( lim*tid, lim + l*tid, n, 0.01 );
                tth[tid] = omp_get_wtime() - tt;

                // We need a few barriers here to ensure that both sides have finish
                // computing, send or retrieving data prior to start the next step
                // First barrier: ensure computations are finished
                #pragma omp barrier

                // First exchanges: retrieving on the host the positions computed on the device
                // while sending to the device the velocities computed on the host
                #pragma omp target update if( tid == 0 ) from( x[0:lim], y[0:lim], z[0:lim] )
                #pragma omp target update if( tid == 1 ) to( vx[lim:l], vy[lim:l], vz[lim:l] )

                // Second barrier: ensure the end of the first two data transfers
                #pragma omp barrier

                // Second exchanges: sending to the device the positions computed on the host
                // while retrieving on the host the velocities computed on the device
                #pragma omp target update if( tid == 0 ) to( x[lim:l], y[lim:l], z[lim:l] )
                #pragma omp target update if( tid == 1 ) from( vx[0:lim], vy[0:lim], vz[0:lim] )
                int a = 0; // useless but here to work around a bug in the Intel 15.0.0 compiler

                // And finally computing the new ratio of work to offload
                #pragma omp single
                ratio = ratio*tth[1] / ((1-ratio)*tth[0] + ratio*tth[1]);
                // There is a implicit final barrier here
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
