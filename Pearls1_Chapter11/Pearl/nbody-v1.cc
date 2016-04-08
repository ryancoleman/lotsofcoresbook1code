// Authors: Gilles Civario and Michael Lysaght
// Copyright 2014 Irish Centre for High-End Computing under MIT license
// See license.txt for the full license of this software
//
// First step of improvement: the code now run on the host CPUs, and offload
// all of its work on the default device (if any). No work is kept on the host.

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

// This extra directives instruct the compiler to also generate the device versions
// of all the variables and functions they enclose
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

void Newton( size_t n, real dt ) {
    const real dtG = dt * G;
    // The processing is done on the default device
    #pragma omp target
    #pragma omp parallel
    {
        #pragma omp for schedule( auto )
        for ( size_t i = 0; i < n; ++i ) {
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
        for ( size_t i = 0; i < n; ++i ) {
            x[i] += vx[i] * dt;
            y[i] += vy[i] * dt;
            z[i] += vz[i] * dt;
        }
    }
}
// End of the "omp declare target" directive
#pragma omp end declare target

int main( int argc, char *argv[] ) {
    const size_t n = 50000;

    bodies( n );

    cout << "Before to start:\n"
         << "  Position of first particle is (" << x[0] << ',' << y[0] << ',' << z[0] << ")\n"
         << "  Position of last particle is (" << x[n-1] << ',' << y[n-1] << ',' << z[n-1] << ")\n";

    double tm = omp_get_wtime();
    // Sending of the necessary input data (to and tofrom) on the device prior to entering
    // the time loop, and retrieving of the output data (tofrom) upon exit of the loop
    #pragma omp target data map( tofrom: x[0:n], y[0:n], z[0:n] ) \
        map( to: vx[0:n], vy[0:n], vz[0:n], m[0:n] )
    for ( int it = 0; it < 100; ++it ) {
        Newton( n, 0.01 );
    }
    tm = omp_get_wtime() - tm;

    cout << "At the end:\n"
         << "  Position of first particle is (" << x[0] << ',' << y[0] << ',' << z[0] << ")\n"
         << "  Position of last particle is (" << x[n-1] << ',' << y[n-1] << ',' << z[n-1] << ")\n";
    cout << "Time was " << tm << "s\n";

    cleanBodies();

    return 0;
}
