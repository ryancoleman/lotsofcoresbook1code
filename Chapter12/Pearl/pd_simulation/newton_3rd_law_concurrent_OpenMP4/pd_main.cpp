//////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014, FLORIAN WENDE, ZIB
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the ZIB nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL FLORIAN WENDE BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
//////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <sched.h>
#include <omp.h>

#include "pd.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////
// DEFINES
/////////////////////////////////////////////////////////////////////////
#define MEASUREMENT
#define BATCH_MODE

#if !defined X
#define X ( 16 )
#endif

#if !defined Y
#define Y ( 16 )
#endif

#if !defined Z
#define Z ( 16 )
#endif

/////////////////////////////////////////////////////////////////////////
// MACROS
/////////////////////////////////////////////////////////////////////////
#define MIN( X, Y ) ( ( X ) < ( Y ) ? ( X ) : ( Y ) )
#define MAX( X, Y ) ( ( X ) > ( Y ) ? ( X ) : ( Y ) )
#define CEIL_N( X, N ) ( ( ( X )/( N )+( ( X )%( N ) ? 1 : 0 ) )*( N ) )

/////////////////////////////////////////////////////////////////////////
// PROTOTYPES
/////////////////////////////////////////////////////////////////////////
static void initParticleChargeMass( const double absMeanCharge, const double meanMass, const uint32_t numParticles, double **charge, double **mass );
static void initParticlePositionVelocity( const double *mass, const double temperature, const uint32_t numParticles_x, const uint32_t numParticles_y, const uint32_t numParticles_z, const double boxExtent, double **position, double **velocity );

/////////////////////////////////////////////////////////////////////////
// MAIN
/////////////////////////////////////////////////////////////////////////
int main( int argc, char **argv ) {

  uint32_t
    numParticles_x = X,
    numParticles_y = Y,
    numParticles_z = Z,
    numParticles = numParticles_x*numParticles_y*numParticles_z;

  double
    absMeanCharge = 1.0,
    meanMass = 1.0,
    epsilon = 1.0,
    sigma = 1.0,
    temperature = 1.0,
    timeStep = 5.0E-5,
    particleDensity = 2.0,
    boxExtent = (double)( (uint32_t)( exp( log( numParticles/particleDensity )/3.0 )+0.5 ) ),
    *position = NULL,
    *velocity = NULL,
    *charge = NULL,
    *mass = NULL;

#if defined USE_HOST
  omp_set_nested( 2 );
#endif

  initParticleChargeMass( absMeanCharge, meanMass, numParticles, &charge, &mass );
  initParticlePositionVelocity( mass, temperature, numParticles_x, numParticles_y, numParticles_z, boxExtent, &position, &velocity );

  pd
    pdSimulation( epsilon, sigma, temperature, boxExtent, numParticles, position, velocity, charge, mass );

  if( position != NULL )
    delete [] position;
  if( velocity != NULL )
    delete [] velocity;
  if( charge != NULL )
    delete [] charge;
  if( mass != NULL )
    delete [] mass;
  
  // the pd-object now can be used for computation
  uint32_t
    numIterations = 500,
    measurements = 10;
  double
    time = omp_get_wtime();

  pdSimulation.computeForce();

  pdSimulation.update( 1.0*timeStep, 0.5*timeStep );

  for( uint32_t i=0; i<numIterations; i++ ) {

    pdSimulation.computeForce();

#if defined MEASUREMENT

    double
      ePot,eKin;

    if( ( i > 0 ) && ( i%measurements ) == 0 ) {

      // same timestep for both position and velocity
      pdSimulation.update( 0.0*timeStep, 0.5*timeStep );
      // compute energy
      pdSimulation.getEnergy( &ePot, &eKin );
      // finalize integration step
      pdSimulation.update( 1.0*timeStep, 0.5*timeStep );

      cout << ::setprecision( 5 ) << ::scientific << ePot << "\t" << eKin << "\t" << ePot+eKin << endl;

    } else {

      pdSimulation.update( 1.0*timeStep, 1.0*timeStep );

    }
    
#else // MEASUREMENT

    pdSimulation.update( 1.0*timeStep, 1.0*timeStep );

#endif // !MEASUREMENT

  } 
  
  time = omp_get_wtime()-time;
  
  // print out information
  cout << ::fixed << "# __INFO: elapsed time " << time*1.0E3 << "ms" << endl;
  cout << "# __INFO: threads used " << pdSimulation.getNumThreads() << endl;
  cout << "# __INFO: gigainteractions " << (uint64_t)( numIterations )*( numParticles*( numParticles-1 ) )*1.0E-9/time << endl;

  // print batch output
#if defined BATCH_MODE
  cout << X*Y*Z << "\t" << pdSimulation.getNumThreads() << "\t" << (uint64_t)( numIterations )*( numParticles*( numParticles-1 ) )*1.0E-9/time << endl;
#endif

  exit( 0 );

}

/////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION
/////////////////////////////////////////////////////////////////////////
static void initParticleChargeMass( const double absMeanCharge, const double meanMass, const uint32_t numParticles, double **charge, double **mass ) {

  // allocate memory
  ( *charge ) = (double *) new double[numParticles];
  ( *mass ) = (double *) new double[numParticles];

  // use drand48() random number generator
  srand48( 1 );

  // mass: randomly distributed around MASS with 10% variation
  cout << "# __INFO: set up particle mass...";

  for( uint32_t i=0; i<numParticles; i++ )
    ( *mass )[i] = meanMass*( 1.0+0.1*( drand48()*2.0-1.0 ) );

  cout << "done" << endl;

  // charge: randomly distributed around (+/-)CHARGE with 10% variation
  cout << "# __INFO: set up particle charge...";

  for( uint32_t i=0; i<numParticles; i++ )
    ( *charge )[i] = absMeanCharge*( 1.0+0.1*( drand48()*2.0-1.0 ) )*( drand48() < 0.5 ? -1.0 : 1.0 );

  cout << "done" << endl;
  
}

static void initParticlePositionVelocity( const double *mass, const double temperature, const uint32_t numParticles_x, const uint32_t numParticles_y, const uint32_t numParticles_z, const double boxExtent, double **position, double **velocity ) {

  uint32_t
    numParticles = numParticles_x*numParticles_y*numParticles_z;

  // allocate memory
  ( *position ) = (double *) new double[3*numParticles];
  ( *velocity ) = (double *) new double[3*numParticles];

  // use drand48() random number generator
  srand48( 2 );

  // positions: regular grid, ( x1,y1,z1 )( x2,y2,z2 )...
  cout << "# __INFO: set up particle position...";

  double
    particleSpacing = (double)( boxExtent )/( MAX( numParticles_x, MAX( numParticles_y, numParticles_z ) )+1 ),
    xShift = ( numParticles_x&1 ? 0.0 : 0.5 );

  for( int32_t z=0; z<numParticles_z; z++ ) {

    double
      zComponent = ( ( -numParticles_z/2 )+z+( numParticles_z&1 ? 0.0 : 0.5 ) )*particleSpacing;
    
    for( int32_t y=0; y<numParticles_y; y++ ) {

      double
	yComponent = ( ( -numParticles_y/2 )+y+( numParticles_y&1 ? 0.0 : 0.5 ) )*particleSpacing;
      
      for( int32_t x=0; x<numParticles_x; x++ ) {

	// x
	( *position )[3*( z*numParticles_y*numParticles_x+y*numParticles_x+x )+0] = ( ( -numParticles_x/2 )+x+xShift )*particleSpacing;
	// y
	( *position )[3*( z*numParticles_y*numParticles_x+y*numParticles_x+x )+1] = yComponent;
	// z
	( *position )[3*( z*numParticles_y*numParticles_x+y*numParticles_x+x )+2] = zComponent;

      }
      
    }

  }

  cout << "done" << endl;

  // velocities: Maxwell distribution with TEMPERATURE, mass[], and with mean velocity zero
  cout << "# __INFO: set up particle velocity...";

  double
    a1,a2,a3,
    b1,b2,b3,
    r1,r2,r3,
    s1,s2,s3,
    scaleFactor,
    meanVelocity[3]; // ( vx,vy,vz )

  meanVelocity[0] = 0.0;
  meanVelocity[1] = 0.0;
  meanVelocity[2] = 0.0;

  for( uint32_t i=0; i<numParticles; i+=2 ) {
    
    do {
      
      a1 = drand48()*2.0-1.0;
      b1 = drand48()*2.0-1.0;
      
      r1 = a1*a1+b1*b1;
      
    } while( r1 >= 1.0 );
    
    do {
      
      a2 = drand48()*2.0-1.0;
      b2 = drand48()*2.0-1.0;
      
      r2 = a2*a2+b2*b2;
      
    } while( r2 >= 1.0 );
    
    do {
      
      a3 = drand48()*2.0-1.0;
      b3 = drand48()*2.0-1.0;
      
      r3 = a3*a3+b3*b3;
      
    } while( r3 >= 1.0 );
    
    s1 = sqrt( -2.0*log( r1 )/r1 );
    s2 = sqrt( -2.0*log( r2 )/r2 );
    s3 = sqrt( -2.0*log( r3 )/r3 );

    scaleFactor = sqrt( temperature/mass[i+0] );

    ( *velocity )[3*i+0] = a1*s1*scaleFactor;
    ( *velocity )[3*i+1] = a2*s2*scaleFactor;
    ( *velocity )[3*i+2] = a3*s3*scaleFactor;
    
    meanVelocity[0] += ( *velocity )[3*i+0];
    meanVelocity[1] += ( *velocity )[3*i+1];
    meanVelocity[2] += ( *velocity )[3*i+2];

    scaleFactor = sqrt( temperature/mass[i+1] );

    ( *velocity )[3*i+3] = b1*s1*scaleFactor;
    ( *velocity )[3*i+4] = b2*s2*scaleFactor;
    ( *velocity )[3*i+5] = b3*s3*scaleFactor;

    meanVelocity[0] += ( *velocity )[3*i+3];
    meanVelocity[1] += ( *velocity )[3*i+4];
    meanVelocity[2] += ( *velocity )[3*i+5];

  }

  meanVelocity[0] /= numParticles;
  meanVelocity[1] /= numParticles;
  meanVelocity[2] /= numParticles;

  for( uint32_t i=0; i<numParticles; i++ ) {

    ( *velocity )[3*i+0] -= meanVelocity[0];
    ( *velocity )[3*i+1] -= meanVelocity[1];
    ( *velocity )[3*i+2] -= meanVelocity[2];

  }

  cout << "done" << endl;

}
