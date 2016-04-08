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
  this->force = (double *) _mm_malloc( NUM_GROUPS*3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ), ALIGNMENT );
  this->charge = (double *) _mm_malloc( CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ), ALIGNMENT );
  this->mass = (double *) _mm_malloc( CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ), ALIGNMENT );

  memset( this->position, 0x0, 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
  memset( this->velocity, 0x0, 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
  memset( this->force, 0x0, NUM_GROUPS*3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
  memset( this->charge, 0x0, CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );
  memset( this->mass, 0x0, CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );

  for( uint32_t i=0; i<NUM_GROUPS; i++ ) {

    this->forceTemp[i] = (double *) _mm_malloc( 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ), ALIGNMENT );
    memset( this->forceTemp[i], 0x0, 3*CEIL_N( numParticles, SIMD_WIDTH )*sizeof( double ) );

  }

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

  // create thread groups on Xeon Phi before hand & one after another
#pragma omp parallel num_threads( NUM_GROUPS )
  {

    uint32_t
      groupId = omp_get_thread_num();

#pragma omp critical ( CREATE_THREAD_GROUPS )
    {

#if !defined USE_HOST
#pragma offload target( mic : device )
      {
#endif

#pragma omp parallel num_threads( NUM_THREADS_PER_GROUP )
	{

#if defined THREAD_PINNING
	  uint32_t
	    ompId = omp_get_thread_num();
	  cpu_set_t
	    cpuMask;

	  CPU_ZERO( &cpuMask );
#if defined USE_CPU
	  CPU_SET( groupId*NUM_THREADS_PER_GROUP+ompId, &cpuMask );
#else
	  CPU_SET( 1+groupId*NUM_THREADS_PER_GROUP+ompId, &cpuMask );
#endif
	  sched_setaffinity( 0, sizeof( cpu_set_t ), &cpuMask );
#else
	  ; // just crete threads
#endif

	}

#if !defined USE_HOST
      }
#endif

    } // critical: CREATE_THREAD_GROUPS

#if !defined USE_HOST
    // transfer per-group data to Xeon Phi
    double
      *_forceTemp = this->forceTemp[groupId];

#pragma offload_transfer target( mic : device ) \
  in( _forceTemp : length( 3*CEIL_N( numParticles, SIMD_WIDTH ) ) align( ALIGNMENT ) alloc_if( true ) free_if( false ) )
#endif

  }

#if !defined USE_HOST
  // offload data: allocation + data transfer
  double
    *_position = this->position,
    *_velocity = this->velocity,
    *_force = this->force,
    *_charge = this->charge,
    *_mass = this->mass;

#pragma offload_transfer target( mic : device ) \
  in( _position : length( 3*CEIL_N( numParticles, SIMD_WIDTH ) ) align( ALIGNMENT ) alloc_if( true ) free_if( false ) ) \
  in( _velocity : length( 3*CEIL_N( numParticles, SIMD_WIDTH ) ) align( ALIGNMENT ) alloc_if( true ) free_if( false ) ) \
  in( _force : length( NUM_GROUPS*3*CEIL_N( numParticles, SIMD_WIDTH ) ) align( ALIGNMENT ) alloc_if( true ) free_if( false ) ) \
  in( _charge : length( CEIL_N( numParticles, SIMD_WIDTH ) ) align( ALIGNMENT ) alloc_if( true ) free_if( false ) ) \
  in( _mass : length( CEIL_N( numParticles, SIMD_WIDTH ) ) align( ALIGNMENT ) alloc_if( true ) free_if( false ) )
#endif

}

pd::~pd() {

#if !defined USE_HOST
  // offload data: release
  double
    *_position = this->position,
    *_velocity = this->velocity,
    *_force = this->force,
    *_charge = this->charge,
    *_mass = this->mass;

#pragma offload target( mic : device ) \
  in( _position : length( 0 ) alloc_if( false ) free_if( true ) ) \
  in( _velocity : length( 0 ) alloc_if( false ) free_if( true ) ) \
  in( _force : length( 0 ) alloc_if( false ) free_if( true ) ) \
  in( _charge : length( 0 ) alloc_if( false ) free_if( true ) ) \
  in( _mass : length( 0 ) alloc_if( false ) free_if( true ) )
  {
    ; // do nothing
  }
#endif

  for( uint32_t i=0; i<NUM_GROUPS; i++ ) {

#if !defined USE_HOST
    double
      *_forceTemp = this->forceTemp[i];

#pragma offload target( mic : device ) \
  in( _forceTemp : length( 0 ) alloc_if( false ) free_if( true ) )
    {
      ; // do nothing
    }
#endif
    
    if( this->forceTemp[i] != NULL )
      _mm_free( this->forceTemp[i] );
    this->forceTemp[i] = NULL;
    
  }

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
