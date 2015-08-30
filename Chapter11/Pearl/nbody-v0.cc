// Authors: Gilles Civario and Michael Lysaght
// Copyright 2014 Irish Centre for High-End Computing under MIT license
// See license.txt for the full license of this software
//
// Initial version of the code, only parallelised and vectorised to run in native
// mode on a host (eg a Xeon processor) or on a device (eg a Xeon Phi co-processor).

#include <iostream>
#include <cstdlib>
#include <cmath>

#include <omp.h>

using namespace std;

// Definition of a generic floating point type, and its corresponding square root function
#ifdef DOUBLE
    typedef double real;
    #define Sqrt( x ) sqrt( x )
#else
    typedef float real;
    #define Sqrt( x ) sqrtf( x )
#endif

// Gravitational constant
const real G = 6.67384e-11;

// Our work arrays for positions (x, y, z), velocities (vx, vy, vz) and masses (m)
real *x, *y, *z, *vx, *vy, *vz, *m;

// Function to allocate and randomly fill an array
void randomFill( real *&x, real inf, real sup, size_t n ) {
    x = new real[n];
    const real mult = (sup - inf) / real( RAND_MAX );
    for ( size_t i = 0; i != n; ++i ) {
        x[i] = rand() * mult + inf;
    }
}

// Function to allocate and fill all our work arrays
void bodies( size_t n ) {
    randomFill( x, -100, 100, n );
    randomFill( y, -100, 100, n );
    randomFill( z, -100, 100, n );
    randomFill( vx, -10, 10, n );
    randomFill( vy, -10, 10, n );
    randomFill( vz, -10, 10, n );
    randomFill( m, 1, 10000, n );    
}

// Function to release all allocated memory
void cleanBodies() {
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] vx;
    delete[] vy;
    delete[] vz;
    delete[] m;
}

// Function computing a time step of Newtonian interactions for our entire system
void Newton( size_t n, real dt ) {
    const real dtG = dt * G;
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

int main( int argc, char *argv[] ) {
    // Number of particles (bodies) in our simulated "universe"
    const size_t n = 50000;

    // Creating the arrays and filling them
    bodies( n );

    // Some printing for checking the validity of the various optimisations
    cout << "Before to start:\n"
         << "  Position of first particle is (" << x[0] << ',' << y[0] << ',' << z[0] << ")\n"
         << "  Position of last particle is (" << x[n-1] << ',' << y[n-1] << ',' << z[n-1] << ")\n";

    double tm = omp_get_wtime();
    // The main loop
    for ( int it = 0; it < 100; ++it ) {
        Newton( n, 0.01 );
    }
    tm = omp_get_wtime() - tm;

    // Final printing of some particle positions. A validity checking is that this shouldn't
    // change with the various optimisations we will do on the code
    cout << "At the end:\n"
         << "  Position of first particle is (" << x[0] << ',' << y[0] << ',' << z[0] << ")\n"
         << "  Position of last particle is (" << x[n-1] << ',' << y[n-1] << ',' << z[n-1] << ")\n";
    // And some timing to see how performance evolves along with the versions of the code
    cout << "Time was " << tm << "s\n";

    // Releasing the memory
    cleanBodies();

    return 0;
}
