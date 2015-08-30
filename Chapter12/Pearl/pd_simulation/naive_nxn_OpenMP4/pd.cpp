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
#include <omp.h>
#include <immintrin.h>

#include "pd.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////
// DEFINES
/////////////////////////////////////////////////////////////////////////
//#define THREAD_PINNING

/////////////////////////////////////////////////////////////////////////
// MACROS
/////////////////////////////////////////////////////////////////////////
#define CEIL_N( X, N ) ( ( ( X )/( N )+( ( X )%( N ) ? 1 : 0 ) )*( N ) )

/////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION
/////////////////////////////////////////////////////////////////////////
pd::pd( const double epsilon,
	const double sigma,
	const double temperature,
	const double boxExtent,
	const uint32_t numParticles,
	const double *position,
	const double *velocity,
	const double *charge,
	const double *mass ) {

  this->epsilon = epsilon;
  this->sigma = sigma;
  this->temperature = temperature;
  this->boxExtent = boxExtent;
  this->numParticles = numParticles;

  // memory allocation:
  this->position = (double *) _mm_malloc( 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ), ALIGNMENT );
  this->velocity = (double *) _mm_malloc( 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ), ALIGNMENT );
  this->force = (double *) _mm_malloc( 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ), ALIGNMENT );
  this->charge = (double *) _mm_malloc( CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ), ALIGNMENT );
  this->mass = (double *) _mm_malloc( CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ), ALIGNMENT );

  memset( this->position, 0x0, 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
  memset( this->velocity, 0x0, 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
  memset( this->force, 0x0, 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
  memset( this->charge, 0x0, CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
  memset( this->mass, 0x0, CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );

  // data format is expected to be (x,y,z)(x,y,z)...
  for( uint32_t i=0; i<numParticles; i++ ) {

    this->position[0*CEIL_N( numParticles, SIMD_WIDTH )+i] = position[3*i+0];
    this->position[1*CEIL_N( numParticles, SIMD_WIDTH )+i] = position[3*i+1];
    this->position[2*CEIL_N( numParticles, SIMD_WIDTH )+i] = position[3*i+2];

    this->velocity[0*CEIL_N( numParticles, SIMD_WIDTH )+i] = velocity[3*i+0];
    this->velocity[1*CEIL_N( numParticles, SIMD_WIDTH )+i] = velocity[3*i+1];
    this->velocity[2*CEIL_N( numParticles, SIMD_WIDTH )+i] = velocity[3*i+2];

    this->charge[i] = charge[i];

    this->mass[i] = mass[i];

  }

#if !defined USE_HOST
  // create thread groups on Xeon Phi before hand
#pragma omp target device( micDevice )
  {
#endif

#pragma omp parallel num_threads( NUM_THREADS )
    {
      
#if defined THREAD_PINNING
      uint32_t
	ompId = omp_get_thread_num();
      cpu_set_t
	cpuMask;
      
      CPU_ZERO( &cpuMask );
#if defined USE_HOST
      CPU_SET( ompId, &cpuMask );
#else
      CPU_SET( 1+ompId, &cpuMask );
#endif
      sched_setaffinity( 0, sizeof( cpu_set_t ), &cpuMask );
#else
      ; // just crete threads
#endif

    }
    
#if !defined USE_HOST
  }
#endif

#if !defined USE_HOST
  // offload data: allocation + data transfer: we use the "pointer-to-uint64_t" trick
  uint32_t
    numParticlesPadded = CEIL_N( numParticles, SIMD_WIDTH );
  uint64_t
    _devPtr_position,
    _devPtr_velocity,
    _devPtr_force,
    _devPtr_charge,
    _devPtr_mass;
  double
    *_position = this->position,
    *_velocity = this->velocity,
    *_force = this->force,
    *_charge = this->charge,
    *_mass = this->mass;

#pragma omp target device( micDevice ) \
  map( to : _position[0:( 3*numParticlesPadded )] ) \
  map( to : _velocity[0:( 3*numParticlesPadded )] ) \
  map( to : _force[0:( 3*numParticlesPadded )] ) \
  map( to : _charge[0:numParticlesPadded] ) \
  map( to : _mass[0:numParticlesPadded] )
  {

    double
      *ptr = NULL;

    // create persistent buffers on the Xeon Phi
    ptr = (double *) _mm_malloc( 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ), ALIGNMENT );
    memcpy( ptr, _position, 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
    _devPtr_position = (uint64_t)ptr;

    ptr = (double *) _mm_malloc( 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ), ALIGNMENT );
    memcpy( ptr, _velocity, 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
    _devPtr_velocity = (uint64_t)ptr;

    ptr = (double *) _mm_malloc( 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ), ALIGNMENT );
    memcpy( ptr, _force, 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
    _devPtr_force = (uint64_t)ptr;

    ptr = (double *) _mm_malloc( CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ), ALIGNMENT );
    memcpy( ptr, _charge, CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
    _devPtr_charge = (uint64_t)ptr;

    ptr = (double *) _mm_malloc( CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ), ALIGNMENT );
    memcpy( ptr, _mass, CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
    _devPtr_mass = (uint64_t)ptr;

    ptr = NULL;

  }

  this->devPtr_position = _devPtr_position;
  this->devPtr_velocity = _devPtr_velocity;
  this->devPtr_force = _devPtr_force;
  this->devPtr_charge = _devPtr_charge;
  this->devPtr_mass = _devPtr_mass;
#endif

}

pd::~pd() {

#if !defined USE_HOST
  // offload data: release
  uint64_t
    _devPtr_position = this->devPtr_position,
    _devPtr_velocity = this->devPtr_velocity,
    _devPtr_force = this->devPtr_force,
    _devPtr_charge = this->devPtr_charge,
    _devPtr_mass = this->devPtr_mass;

#pragma omp target device( micDevice )
  {
    _mm_free( (double *)_devPtr_position );
    _mm_free( (double *)_devPtr_velocity );
    _mm_free( (double *)_devPtr_force );
    _mm_free( (double *)_devPtr_charge );
    _mm_free( (double *)_devPtr_mass );
  }
#endif

  // release host memory
  if( position != NULL )
    _mm_free( position );
  if( velocity != NULL )
    _mm_free( velocity );
  if( force != NULL )
    _mm_free( force );
  if( charge != NULL )
    _mm_free( charge );
  if( mass != NULL )
    _mm_free( mass );

  position = NULL;
  velocity = NULL;
  force = NULL;
  charge = NULL;
  mass = NULL;

}
