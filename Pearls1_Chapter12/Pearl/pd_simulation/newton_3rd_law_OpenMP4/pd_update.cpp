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
#include <omp.h>

#include "pd.h"

/////////////////////////////////////////////////////////////////////////
//
// update positions and velocities: leap frog integrator
//
//     x[i](n+1) = x[i](n)+v[i](n+1/2)*t_x
//     v[i](n+3/2) = v[i](n+1/2)+( F[i](n+1)/m[i] )*t_v
//
// note: first update velocity then position
//
//       compute F(x(0))
//     ||v(0)   -> v(1/2)        : t_x=1.0, t_v=0.5
//     |+----------------
//     ||x(0)   -> x(1)
//       compute F(x(1))
//     ||v(1/2) -> v(3/2)        : t_x=1.0, t_v=1.0
//     |+----------------
//     ||x(1)   -> x(2)
//       compute F(x(2))
//     ||v(3/2) -> v(5/2)        : t_x=1.0, t_v=1.0
//     |+----------------
//     ||x(2)   -> x(3)
//       ...
//
///////////////////////////////////////////////////////////////////////// 

/////////////////////////////////////////////////////////////////////////
// MACROS
/////////////////////////////////////////////////////////////////////////
#define MIN( X, Y ) ( ( X ) < ( Y ) ? ( X ) : ( Y ) )
#define CEIL_N( X, N ) ( ( ( X )/( N )+( ( X )%( N ) ? 1 : 0 ) )*( N ) )

/////////////////////////////////////////////////////////////////////////
// PROTOTYPES
/////////////////////////////////////////////////////////////////////////
__attribute__((target( mic ))) 
void _update( double *position,
	      double *velocity,
	      const double *force,
	      const double *mass,
	      const double t_x,
	      const double t_v,
	      const double boxExtent,
	      const uint32_t numParticles );

/////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION
///////////////////////////////////////////////////////////////////////// 
void pd::update( const double t_x, const double t_v ) {

  uint32_t
    _numParticles = this->numParticles;
  uint64_t
    _devPtr_position = this->devPtr_position,
    _devPtr_velocity = this->devPtr_velocity,
    _devPtr_force = this->devPtr_force,
    _devPtr_mass = this->devPtr_mass;
  double
    _boxExtent = this->boxExtent;

#if !defined USE_HOST
#pragma omp target device( micDevice )
  {
#endif

    double
      *_position = (double *)_devPtr_position,
      *_velocity = (double *)_devPtr_velocity,
      *_force = (double *)_devPtr_force,
      *_mass = (double *)_devPtr_mass;

    _update( _position, _velocity, _force, _mass, t_x, t_v, _boxExtent, _numParticles );

#if !defined USE_HOST
  }
#endif

}

__attribute__((target( mic ))) 
void _update( double *position,
	      double *velocity,
	      const double *force,
	      const double *mass,
	      const double t_x,
	      const double t_v,
	      const double boxExtent,
	      const uint32_t numParticles ) {

  __assume_aligned( force, ALIGNMENT );
  __assume_aligned( mass, ALIGNMENT );
  __assume_aligned( position, ALIGNMENT );
  __assume_aligned( velocity, ALIGNMENT );
  
  double
    plusBoxExtent = boxExtent,
    minusBoxExtent = -boxExtent,
    plusHalfBoxExtent = 0.5*boxExtent,
    minusHalfBoxExtent = -0.5*boxExtent;

#pragma omp parallel num_threads( NUM_THREADS )
  {

    uint32_t
      ompId = omp_get_thread_num(),
      start = ompId*CEIL_N( CEIL_N( numParticles, NUM_THREADS )/NUM_THREADS, SIMD_WIDTH ),
      stop = MIN( numParticles, start+CEIL_N( CEIL_N( numParticles, NUM_THREADS )/NUM_THREADS, SIMD_WIDTH ) );
    double
      *pos_x = (double *)&position[0*CEIL_N( numParticles, SIMD_WIDTH )],
      *pos_y = (double *)&position[1*CEIL_N( numParticles, SIMD_WIDTH )],
      *pos_z = (double *)&position[2*CEIL_N( numParticles, SIMD_WIDTH )],
      *vel_x = (double *)&velocity[0*CEIL_N( numParticles, SIMD_WIDTH )],
      *vel_y = (double *)&velocity[1*CEIL_N( numParticles, SIMD_WIDTH )],
      *vel_z = (double *)&velocity[2*CEIL_N( numParticles, SIMD_WIDTH )],
      *frc_x = (double *)&force[0*CEIL_N( numParticles, SIMD_WIDTH )],
      *frc_y = (double *)&force[1*CEIL_N( numParticles, SIMD_WIDTH )],
      *frc_z = (double *)&force[2*CEIL_N( numParticles, SIMD_WIDTH )];

    __assume_aligned( pos_x, ALIGNMENT );
    __assume_aligned( pos_y, ALIGNMENT );
    __assume_aligned( pos_z, ALIGNMENT );

    __assume_aligned( vel_x, ALIGNMENT );
    __assume_aligned( vel_y, ALIGNMENT );
    __assume_aligned( vel_z, ALIGNMENT );

    __assume_aligned( frc_x, ALIGNMENT );
    __assume_aligned( frc_y, ALIGNMENT );
    __assume_aligned( frc_z, ALIGNMENT );
    
#pragma vector aligned novecremainder
    for( uint32_t i=start; i<stop; i++ ) {
      
      vel_x[i] += ( frc_x[i]/mass[i] )*t_v;
      vel_y[i] += ( frc_y[i]/mass[i] )*t_v;
      vel_z[i] += ( frc_z[i]/mass[i] )*t_v;
      
    }
    
#pragma vector aligned novecremainder
    for( uint32_t i=start; i<stop; i++ ) {
      
      pos_x[i] += vel_x[i]*t_x;
      pos_y[i] += vel_y[i]*t_x;
      pos_z[i] += vel_z[i]*t_x;
      
    }
   
#pragma vector aligned novecremainder
    for( uint32_t i=start; i<stop; i++ ) { 
      
      if( pos_x[i] < minusHalfBoxExtent ) {
	
	pos_x[i] -= ( plusBoxExtent+2.0*pos_x[i] );
	vel_x[i] *= -1.0;
	
      }
      
      if( pos_x[i] > plusHalfBoxExtent ) {
	
	pos_x[i] -= ( minusBoxExtent+2.0*pos_x[i] );
	vel_x[i] *= -1.0;
	
      }
      
      if( pos_y[i] < minusHalfBoxExtent ) {
	
	pos_y[i] -= ( plusBoxExtent+2.0*pos_y[i] );
	vel_y[i] *= -1.0;
	
      }
      
      if( pos_y[i] > plusHalfBoxExtent ) {
	
      pos_y[i] -= ( minusBoxExtent+2.0*pos_y[i] );
      vel_y[i] *= -1.0;
      
      }
      
      if( pos_z[i] < minusHalfBoxExtent ) {
	
	pos_z[i] -= ( plusBoxExtent+2.0*pos_z[i] );
	vel_z[i] *= -1.0;
	
      }
      
      if( pos_z[i] > plusHalfBoxExtent ) {
	
	pos_z[i] -= ( minusBoxExtent+2.0*pos_z[i] );
	vel_z[i] *= -1.0;
	
      }

    }
    
  }
  
}
