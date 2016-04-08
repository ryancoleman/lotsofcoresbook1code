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
		    double *forceTemp,
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
  double
    *_position = this->position,
    *_force = this->force,
    *_forceTemp = this->forceTemp,
    *_charge = this->charge,
    _sigma = this->sigma,
    _epsilon = this->epsilon;

#if !defined USE_HOST
#pragma offload target( mic : device ) \
  in( _position : length( 0 ) alloc_if( false ) free_if( false ) ) \
  in( _force : length( 0 ) alloc_if( false ) free_if( false ) ) \
  in( _forceTemp : length( 0 ) alloc_if( false ) free_if( false ) ) \
  in( _charge : length( 0 ) alloc_if( false ) free_if( false ) )
  {
#endif

    _computeForce( _position, _force, _forceTemp, _charge, _sigma, _epsilon, _numParticles );

#if !defined USE_HOST
  }
#endif

}

__attribute__((target( mic ))) 
void _computeForce( const double *position,
		    double *force,
		    double *forceTemp,
		    const double *charge,
		    const double sigma,
		    const double epsilon,
		    const uint32_t numParticles ) {

  double
    sigma_6 = sigma*sigma*sigma*sigma*sigma*sigma,
    epsilon_sigma_6 = epsilon*sigma_6,
    inverse_4_pi_epsilon = 1.0,
    *pos_x = (double *)&position[0*CEIL_N( numParticles, SIMD_WIDTH )],
    *pos_y = (double *)&position[1*CEIL_N( numParticles, SIMD_WIDTH )],
    *pos_z = (double *)&position[2*CEIL_N( numParticles, SIMD_WIDTH )],
    *frc_x = (double *)&force[0*CEIL_N( numParticles, SIMD_WIDTH )],
    *frc_y = (double *)&force[1*CEIL_N( numParticles, SIMD_WIDTH )],
    *frc_z = (double *)&force[2*CEIL_N( numParticles, SIMD_WIDTH )],
    *frc_j_x = (double *)&forceTemp[0*CEIL_N( numParticles, SIMD_WIDTH )],
    *frc_j_y = (double *)&forceTemp[1*CEIL_N( numParticles, SIMD_WIDTH )],
    *frc_j_z = (double *)&forceTemp[2*CEIL_N( numParticles, SIMD_WIDTH )];
  
  __assume_aligned( position, ALIGNMENT );
  __assume_aligned( charge, ALIGNMENT );
  __assume_aligned( force, ALIGNMENT );
  __assume_aligned( forceTemp, ALIGNMENT );
  
  __assume_aligned( pos_x, ALIGNMENT );
  __assume_aligned( pos_y, ALIGNMENT );
  __assume_aligned( pos_z, ALIGNMENT );
  
  __assume_aligned( frc_x, ALIGNMENT );
  __assume_aligned( frc_y, ALIGNMENT );
  __assume_aligned( frc_z, ALIGNMENT );
  
  __assume_aligned( frc_j_x, ALIGNMENT );
  __assume_aligned( frc_j_y, ALIGNMENT );
  __assume_aligned( frc_j_z, ALIGNMENT );
 
  memset( frc_x, 0x0, CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
  memset( frc_y, 0x0, CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
  memset( frc_z, 0x0, CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
  memset( frc_j_x, 0x0, CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
  memset( frc_j_y, 0x0, CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
  memset( frc_j_z, 0x0, CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
  
#pragma omp parallel num_threads( NUM_THREADS )
  {

    uint32_t
      ompId = omp_get_thread_num(),
      j;
    double
      r_ij_x,r_ij_y,r_ij_z,
      r_ij_1,r_ij_2,r_ij_3,r_ij_6,r_ij_8,
      f_i_x[UNROLL_FACTOR_Y] __attribute__((align( ALIGNMENT ))),
      f_i_y[UNROLL_FACTOR_Y] __attribute__((align( ALIGNMENT ))),
      f_i_z[UNROLL_FACTOR_Y] __attribute__((align( ALIGNMENT ))),
      f_j_x_temp[UNROLL_FACTOR_X] __attribute__((align( ALIGNMENT ))),
      f_j_y_temp[UNROLL_FACTOR_X] __attribute__((align( ALIGNMENT ))),
      f_j_z_temp[UNROLL_FACTOR_X] __attribute__((align( ALIGNMENT ))),
      f;

    for( uint32_t k=0; k<NUM_THREADS; k++ ) {
    
      for( uint32_t i=( ompId*UNROLL_FACTOR_Y ); i<numParticles; i+=( NUM_THREADS*UNROLL_FACTOR_Y ) ) {

	for( uint32_t ii=0; ii<UNROLL_FACTOR_Y; ii++ ) {

	  f_i_x[ii] = 0.0;
	  f_i_y[ii] = 0.0;
	  f_i_z[ii] = 0.0;

	}

	uint32_t
	  jStart = ( ( k+ompId )%NUM_THREADS )*UNROLL_FACTOR_X,
	  jStop = ( i/UNROLL_FACTOR_X )*UNROLL_FACTOR_X,
	  jIncrement = NUM_THREADS*UNROLL_FACTOR_X;

	for( j=jStart; j<jStop; j+=jIncrement ) {

#pragma vector aligned novecremainder
	  for( uint32_t jj=0; jj<UNROLL_FACTOR_X; jj++ ) {

	    f_j_x_temp[jj] = 0.0;
	    f_j_y_temp[jj] = 0.0;
	    f_j_z_temp[jj] = 0.0;

	  }

	  for( uint32_t ii=0; ii<UNROLL_FACTOR_Y; ii++ ) {

#pragma vector aligned novecremainder
	    for( uint32_t jj=0; jj<UNROLL_FACTOR_X; jj++ ) {
	     
	      r_ij_x = pos_x[i+ii]-pos_x[j+jj];
	      r_ij_y = pos_y[i+ii]-pos_y[j+jj];
	      r_ij_z = pos_z[i+ii]-pos_z[j+jj];
	      
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
	      
	      f = ( 24.0*epsilon_sigma_6*r_ij_8 )*( 2.0*sigma_6*r_ij_6-1.0 )+( inverse_4_pi_epsilon*charge[i+ii]*charge[j+jj]*r_ij_3 );
	      
	      f_i_x[ii] += f*( r_ij_x );
	      f_i_y[ii] += f*( r_ij_y );
	      f_i_z[ii] += f*( r_ij_z );

	      f_j_x_temp[jj] += f*( r_ij_x );
	      f_j_y_temp[jj] += f*( r_ij_y );
	      f_j_z_temp[jj] += f*( r_ij_z );
	      
	    }

	  }

#pragma vector aligned novecremainder 
	  for( uint32_t jj=0; jj<UNROLL_FACTOR_X; jj++ ) {

	    frc_j_x[j+jj] -= f_j_x_temp[jj];
	    frc_j_y[j+jj] -= f_j_y_temp[jj];
	    frc_j_z[j+jj] -= f_j_z_temp[jj];

	  }
	  
	}

	jStart = j;
	
	for( uint32_t ii=0; ii<UNROLL_FACTOR_Y; ii++ ) {

#pragma vector aligned vecremainder
	  for( j=jStart; j<( i+ii ); j++ ) {
	  
	    r_ij_x = pos_x[i+ii]-pos_x[j];
	    r_ij_y = pos_y[i+ii]-pos_y[j];
	    r_ij_z = pos_z[i+ii]-pos_z[j];
	    
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
	    
	    f = ( 24.0*epsilon_sigma_6*r_ij_8 )*( 2.0*sigma_6*r_ij_6-1.0 )+( inverse_4_pi_epsilon*charge[i+ii]*charge[j]*r_ij_3 );
	    
	    f_i_x[ii] += f*( r_ij_x );
	    f_i_y[ii] += f*( r_ij_y );
	    f_i_z[ii] += f*( r_ij_z );
	    
	    frc_j_x[j] -= f*( r_ij_x );
	    frc_j_y[j] -= f*( r_ij_y );
	    frc_j_z[j] -= f*( r_ij_z );

	  }
	  
	}

#pragma vector aligned novecremainder
	for( uint32_t ii=0; ii<UNROLL_FACTOR_Y; ii++ ) {

	  frc_x[i+ii] += f_i_x[ii];
	  frc_y[i+ii] += f_i_y[ii];
	  frc_z[i+ii] += f_i_z[ii];

	}
	
      }

#pragma omp barrier

    }

    uint32_t
      start = ompId*CEIL_N( numParticles, NUM_THREADS*SIMD_WIDTH )/NUM_THREADS,
      stop = MIN( numParticles, start+CEIL_N( numParticles, NUM_THREADS*SIMD_WIDTH )/NUM_THREADS );

#pragma vector aligned novecremainder
    for( uint32_t i=start; i<stop; i++ ) 
      frc_x[i] += frc_j_x[i];

#pragma vector aligned novecremainder
    for( uint32_t i=start; i<stop; i++ )
      frc_y[i] += frc_j_y[i];

#pragma vector aligned novecremainder
    for( uint32_t i=start; i<stop; i++ )
      frc_z[i] += frc_j_z[i];

  }

}
