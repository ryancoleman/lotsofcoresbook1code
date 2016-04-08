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

#include <cstdint>
#include <cstring>
#include <cmath>
#include <omp.h>

#include "pd.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////
//
// compute force: Lennard-Jones + Coulomb
//
//     F[i] = SUM_j { ( 24.0 * epsilon * sigma^6 / r^8 )*( 2.0 * sigma^6 / r^6 - 1.0 )
//                   +( 1.0 / ( 4.0 * pi * epsilon_0 ) )*( q[i] * q[j] / r^3 ) } * ( r[i]-r[j] )
//
//     with: r = || r[i]-r[j] ||
//
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
// DEFINES
/////////////////////////////////////////////////////////////////////////
#define UNROLL_FACTOR_X ( SIMD_WIDTH )
#define UNROLL_FACTOR_Y ( SIMD_WIDTH )

/////////////////////////////////////////////////////////////////////////
// MACROS
/////////////////////////////////////////////////////////////////////////
#define MIN( X, Y ) ( ( X ) < ( Y ) ? ( X ) : ( Y ) )
#define MAX( X, Y ) ( ( X ) > ( Y ) ? ( X ) : ( Y ) )
#define CEIL_N( X, N ) ( ( ( X )/( N )+( ( X )%( N ) ? 1 : 0 ) )*( N ) )
#define FLOOR_N( X, N ) ( ( ( X )/( N ) )*( N ) )

/////////////////////////////////////////////////////////////////////////
// PROTOTYPES
/////////////////////////////////////////////////////////////////////////
__attribute__((target( mic ))) 
void _computeForce( const double *position,
		    double *force,
		    const double *charge,
		    const double sigma,
		    const double epsilon,
		    const uint32_t numParticles );

/////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION
/////////////////////////////////////////////////////////////////////////
void pd::computeForce() {

  uint32_t
    _numParticles = this->numParticles;
  uint64_t
    _devPtr_position = this->devPtr_position,
    _devPtr_force = this->devPtr_force,
    _devPtr_charge = this->devPtr_charge;
  double
    _sigma = this->sigma,
    _epsilon = this->epsilon;

#if !defined USE_HOST
#pragma omp target device( micDevice )
  {
#endif

    double
      *_position = (double *)_devPtr_position,
      *_force = (double *)_devPtr_force,
      *_charge = (double *)_devPtr_charge;

    _computeForce( _position, _force, _charge, _sigma, _epsilon, _numParticles );

#if !defined USE_HOST
  }
#endif

}

__attribute__((target( mic ))) 
void _computeForce( const double *position,
		    double *force,
		    const double *charge,
		    const double sigma,
		    const double epsilon,
		    const uint32_t numParticles ) {

  __assume_aligned( position, ALIGNMENT );
  __assume_aligned( charge, ALIGNMENT );
  __assume_aligned( force, ALIGNMENT );

  double
    sigma_6 = sigma*sigma*sigma*sigma*sigma*sigma,
    epsilon_sigma_6 = epsilon*sigma_6,
    inverse_4_pi_epsilon = 1.0;

#pragma omp parallel num_threads( NUM_THREADS )
  {

    double
      r_ij_x,r_ij_y,r_ij_z,
      r_ij_1,r_ij_2,r_ij_3,r_ij_6,r_ij_8,
      f_i_x,f_i_y,f_i_z,
      f,
      *pos_x = (double *)&position[0*CEIL_N( numParticles, SIMD_WIDTH )],
      *pos_y = (double *)&position[1*CEIL_N( numParticles, SIMD_WIDTH )],
      *pos_z = (double *)&position[2*CEIL_N( numParticles, SIMD_WIDTH )],
      *frc_x = (double *)&force[0*CEIL_N( numParticles, SIMD_WIDTH )],
      *frc_y = (double *)&force[1*CEIL_N( numParticles, SIMD_WIDTH )],
      *frc_z = (double *)&force[2*CEIL_N( numParticles, SIMD_WIDTH )];

#pragma omp for
    for( uint32_t i=0; i<numParticles; i++ ) {

      f_i_x = 0.0;
      f_i_y = 0.0;
      f_i_z = 0.0;

#if defined __AVX__
#pragma vector aligned
#else
#pragma simd
#endif
      for( uint32_t j=0; j<i; j++ ) {

	r_ij_x = pos_x[i]-pos_x[j];
	r_ij_y = pos_y[i]-pos_y[j];
	r_ij_z = pos_z[i]-pos_z[j];

#if defined FAST_MATH
	r_ij_1 = 1.0/sqrtf( r_ij_x*r_ij_x+r_ij_y*r_ij_y+r_ij_z*r_ij_z );
#elif defined FAST_FAST_MATH
	r_ij_1 = (double)( 1.0F/sqrtf( (float)( r_ij_x*r_ij_x+r_ij_y*r_ij_y+r_ij_z*r_ij_z ) ) );
#else
	r_ij_1 = 1.0/sqrt( r_ij_x*r_ij_x+r_ij_y*r_ij_y+r_ij_z*r_ij_z );
#endif		

	r_ij_2 = r_ij_1*r_ij_1;
	r_ij_3 = r_ij_1*r_ij_2;
	r_ij_6 = r_ij_3*r_ij_3;
	r_ij_8 = r_ij_2*r_ij_6;

	f = ( 24.0*epsilon_sigma_6*r_ij_8 )*( 2.0*sigma_6*r_ij_6-1.0 )+( inverse_4_pi_epsilon*charge[i]*charge[j]*r_ij_3 );
	
	f_i_x += f*( r_ij_x );
	f_i_y += f*( r_ij_y );
	f_i_z += f*( r_ij_z );
	
      }

#if defined __AVX__
#pragma vector aligned
#else
#pragma simd
#endif
      for( uint32_t j=( i+1 ); j<numParticles; j++ ) {

	r_ij_x = pos_x[i]-pos_x[j];
	r_ij_y = pos_y[i]-pos_y[j];
	r_ij_z = pos_z[i]-pos_z[j];

#if defined FAST_MATH
	r_ij_1 = 1.0/sqrtf( r_ij_x*r_ij_x+r_ij_y*r_ij_y+r_ij_z*r_ij_z );
#elif defined FAST_FAST_MATH
	r_ij_1 = (double)( 1.0F/sqrtf( (float)( r_ij_x*r_ij_x+r_ij_y*r_ij_y+r_ij_z*r_ij_z ) ) );
#else
	r_ij_1 = 1.0/sqrt( r_ij_x*r_ij_x+r_ij_y*r_ij_y+r_ij_z*r_ij_z );
#endif		

	r_ij_2 = r_ij_1*r_ij_1;
	r_ij_3 = r_ij_1*r_ij_2;
	r_ij_6 = r_ij_3*r_ij_3;
	r_ij_8 = r_ij_2*r_ij_6;

	f = ( 24.0*epsilon_sigma_6*r_ij_8 )*( 2.0*sigma_6*r_ij_6-1.0 )+( inverse_4_pi_epsilon*charge[i]*charge[j]*r_ij_3 );
	
	f_i_x += f*( r_ij_x );
	f_i_y += f*( r_ij_y );
	f_i_z += f*( r_ij_z );
	
      }

      frc_x[i] = f_i_x;
      frc_y[i] = f_i_y;
      frc_z[i] = f_i_z;

    }

  }

}
