#//////////////////////////////////////////////////////////////////////////////////
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

#if !defined PD_H
#define PD_H

#include <cstdint>
#include <sched.h>

#define ALIGNMENT ( 64 )
#define SIMD_WIDTH ( 8 )

//#define USE_HOST

#if !defined NUM_THREADS
#if defined USE_HOST
#define NUM_THREADS ( 32 )
#else
#define NUM_THREADS ( 240 )
#endif
#endif

class pd {

 public:
  
  pd( const double epsilon,
      const double sigma,
      const double temperature,
      const double boxExtent,
      const uint32_t numParticles,
      const double *position,
      const double *velocity,
      const double *charge,
      const double *mass );
  
  ~pd();

  void computeForce();

  void update( const double t_x, const double t_v );

  void getEnergy( double *ePot, double *eKin );

  void getPosition( double *position ) const { ; /*** not implemented yet ***/ };
  
  void getVelocity( double *velocity ) const { ; /*** not implemented yet ***/ };

  uint32_t getNumThreads() const { return NUM_THREADS; }
  
 private:
  
  uint32_t
    numParticles,
    micDevice = 0;
  double 
    epsilon,
    sigma,
    temperature,
    boxExtent,
    *position = NULL,
    *velocity = NULL,
    *force = NULL,
    *charge = NULL,
    *mass = NULL;
  uint64_t
    devPtr_position,
    devPtr_velocity,
    devPtr_force,
    devPtr_charge,
    devPtr_mass;

};

#endif
