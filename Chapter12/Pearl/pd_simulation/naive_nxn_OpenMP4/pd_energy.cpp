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
#include <cmath>
#include <omp.h>

#include "pd.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////
//
// compute potential and kinetic energy:
//
//     E_pot[i] = SUM_j { ( 4.0 * epsilon * sigma^6 / r^6 )*( sigma^6 / r^6 - 1.0 ) 
//                       +( 1.0 / ( 4.0 * pi * epsilon_0 ) )*( q[i] * q[j] / r ) }
//
//     E_kin[i] = m[i] * v[i]^2 / 2.0
//
//     with: r = || r[i]-r[j] ||
//
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
// MACROS
/////////////////////////////////////////////////////////////////////////
#define CEIL_N( X, N ) ( ( ( X )/( N )+( ( X )%( N ) ? 1 : 0 ) )*( N ) )

/////////////////////////////////////////////////////////////////////////
// DATATYPES
/////////////////////////////////////////////////////////////////////////
__attribute__((target( mic )))
typedef struct {

  double ePot;
  double eKin;

} energy_t;

/////////////////////////////////////////////////////////////////////////
// PROTOTYPES
/////////////////////////////////////////////////////////////////////////
__attribute__((target( mic ))) 
energy_t _getEnergy( const double *position,
		     const double *velocity,
		     const double *charge,
		     const double *mass,
		     const double sigma,
		     const double epsilon,
		     const uint32_t numParticles );

/////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION
/////////////////////////////////////////////////////////////////////////
void pd::getEnergy( double *ePot, double *eKin ) {

  uint32_t
    _numParticles = this->numParticles;
  uint64_t
    _devPtr_position = this->devPtr_position,
    _devPtr_velocity = this->devPtr_velocity,
    _devPtr_charge = this->devPtr_charge,
    _devPtr_mass = this->devPtr_mass;
  double
    _sigma = this->sigma,
    _epsilon = this->epsilon,
    _ePot,
    _eKin;

#if !defined USE_HOST
#pragma omp target device( micDevice )
  {
#endif

    double
      *_position = (double *)_devPtr_position,
      *_velocity = (double *)_devPtr_velocity,
      *_charge = (double *)_devPtr_charge,
      *_mass = (double *)_devPtr_mass;
    energy_t
      energy;

    energy = _getEnergy( _position, _velocity, _charge, _mass, _sigma, _epsilon, _numParticles );

    _ePot = energy.ePot;
    _eKin = energy.eKin;
   
#if !defined USE_HOST
  }
#endif

  ( *ePot ) = _ePot;
  ( *eKin ) = _eKin;

}

__attribute__((target( mic ))) 
energy_t _getEnergy( const double *position,
		     const double *velocity,
		     const double *charge,
		     const double *mass,
		     const double sigma,
		     const double epsilon,
		     const uint32_t numParticles ) {
  
  __assume_aligned( position, ALIGNMENT );
  __assume_aligned( velocity, ALIGNMENT );
  __assume_aligned( charge, ALIGNMENT );
  __assume_aligned( mass, ALIGNMENT );

  double
    sigma_6 = sigma*sigma*sigma*sigma*sigma*sigma,
    epsilon_sigma_6 = epsilon*sigma_6,
    inverse_4_pi_epsilon = 1.0,
    _ePot,
    _eKin;
  energy_t
    energy;

  _ePot = 0.0;
  _eKin = 0.0;

#pragma omp parallel num_threads( NUM_THREADS ) reduction( + : _ePot ) reduction( + : _eKin )
  {

    double
      r_ij_x,r_ij_y,r_ij_z,
      r_ij_1,r_ij_2,r_ij_3,r_ij_6,
      *pos_x = (double *)&position[0*CEIL_N( numParticles, SIMD_WIDTH )],
      *pos_y = (double *)&position[1*CEIL_N( numParticles, SIMD_WIDTH )],
      *pos_z = (double *)&position[2*CEIL_N( numParticles, SIMD_WIDTH )],
      *vel_x = (double *)&velocity[0*CEIL_N( numParticles, SIMD_WIDTH )],
      *vel_y = (double *)&velocity[1*CEIL_N( numParticles, SIMD_WIDTH )],
      *vel_z = (double *)&velocity[2*CEIL_N( numParticles, SIMD_WIDTH )];

    __assume_aligned( pos_x, ALIGNMENT );
    __assume_aligned( pos_y, ALIGNMENT );
    __assume_aligned( pos_z, ALIGNMENT );
    __assume_aligned( vel_x, ALIGNMENT );
    __assume_aligned( vel_y, ALIGNMENT );
    __assume_aligned( vel_z, ALIGNMENT );


#pragma omp for
    for( uint32_t i=0; i<numParticles; i++ ) {

#pragma vector aligned
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
	r_ij_3 = r_ij_1*r_ij_1*r_ij_1;
	r_ij_6 = r_ij_3*r_ij_3;

	_ePot += ( 4.0*epsilon_sigma_6*r_ij_6 )*( sigma_6*r_ij_6-1.0 )+( inverse_4_pi_epsilon*charge[i]*charge[j]*r_ij_1 );
	
      }

      _eKin += 0.5*mass[i]*( vel_x[i]*vel_x[i]+vel_y[i]*vel_y[i]+vel_z[i]*vel_z[i] );
      
    }
    
  }

  energy.ePot = _ePot;
  energy.eKin = _eKin;

  return energy;

}
