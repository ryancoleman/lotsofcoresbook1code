/***************************************************************************
 * Copyright (2012)2 (03-2014)3 Intel Corporation All Rights Reserved.
 *
 * The source code contained or described herein and all documents related to 
 * the source code ("Material") are owned by Intel Corporation or its suppliers 
 * or licensors. Title to the Material remains with Intel Corporation or its 
 * suppliers and licensors. The Material contains trade secrets and proprietary 
 * and confidential information of Intel or its suppliers and licensors. The 
 * Material is protected by worldwide copyright and trade secret laws and 
 * treaty provisions. No part of the Material may be used, copied, reproduced, 
 * modified, published, uploaded, posted, transmitted, distributed, or disclosed 
 * in any way without Intel’s prior express written permission.
 *
 * No license under any patent, copyright, trade secret or other intellectual 
 * property right is granted to or conferred upon you by disclosure or delivery 
 * of the Materials, either expressly, by implication, inducement, estoppel or 
 * otherwise. Any license under such intellectual property rights must be express 
 * and approved by Intel in writing.
 * ***************************************************************************/

/*****************************************************************************
 * ! Content:
 * ! Implementation example of ISO-3DFD implementation for 
 * !   Intel(R) Xeon Phi(TM) and Intel(R) Xeon.
 * ! Version 07
 * ! leonardo.borges@intel.com
 * ! cedric.andreolli@intel.com
 * !****************************************************************************/

#ifndef _TOOLS_INCLUDE
#define _TOOLS_INCLUDE

#include <stddef.h>
#include <sys/time.h>
#include <time.h>
#include "iso-3dfd.h"

// NOTE: the use of clock_gettime() below requires you to link
// with -lrt (the real-time clock)
double walltime() // seconds
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);

    return ((double)(ts.tv_sec) + 
            1e-09 * (double)(ts.tv_nsec));
}

void first_touch(float* tab_base, long n1, long n2, long n3, long n1_Tblock, long n2_Tblock, long n3_Tblock, long num_threads){
	long n3End = n3 - HALF_LENGTH;
	long n2End = n2 - HALF_LENGTH;
	long n1End = n1 - HALF_LENGTH;
	long dimn1n2 = n1*n2;
	#pragma omp parallel OMP_N_THREADS
	{
		float* tab;
	#ifdef BLOCK_X_Y_Z
		#pragma omp for OMP_SCHEDULE collapse(3)	
		for(long bx=HALF_LENGTH; bx<n1End; bx+=n1_Tblock){	
			for(long by=HALF_LENGTH; by<n2End; by+=n2_Tblock){	
				for(long bz=HALF_LENGTH; bz<n3End; bz+=n3_Tblock){	
	#endif

	#ifdef BLOCK_X_Z_Y
		#pragma omp for OMP_SCHEDULE collapse(3)	
		for(long bx=HALF_LENGTH; bx<n1End; bx+=n1_Tblock){	
			for(long bz=HALF_LENGTH; bz<n3End; bz+=n3_Tblock){	
				for(long by=HALF_LENGTH; by<n2End; by+=n2_Tblock){	
	#endif

	#ifdef BLOCK_Y_Z_X
	        #pragma omp for OMP_SCHEDULE collapse(3) 
	        for(long by=HALF_LENGTH; by<n2End; by+=n2_Tblock){
	                for(long bz=HALF_LENGTH; bz<n3End; bz+=n3_Tblock){       
	        		for(long bx=HALF_LENGTH; bx<n1End; bx+=n1_Tblock){
	#endif

	#ifdef BLOCK_Y_X_Z
	        #pragma omp for OMP_SCHEDULE collapse(3) 
	        for(long by=HALF_LENGTH; by<n2End; by+=n2_Tblock){
	                for(long bx=HALF_LENGTH; bx<n1End; bx+=n1_Tblock){
	                	for(long bz=HALF_LENGTH; bz<n3End; bz+=n3_Tblock){    
	#endif  

	#ifdef BLOCK_Z_X_Y
	        #pragma omp for OMP_SCHEDULE collapse(3) 
	        for(long bz=HALF_LENGTH; bz<n3End; bz+=n3_Tblock){
	                for(long bx=HALF_LENGTH; bx<n1End; bx+=n1_Tblock){
	        		for(long by=HALF_LENGTH; by<n2End; by+=n2_Tblock){
	#endif 

	#ifdef BLOCK_Z_Y_X
	        #pragma omp for collapse(3) 
	       	for(long bz=HALF_LENGTH; bz<n3End; bz+=n3_Tblock){
	        	for(long by=HALF_LENGTH; by<n2End; by+=n2_Tblock){
	                	for(long bx=HALF_LENGTH; bx<n1End; bx+=n1_Tblock){
	#endif 
//											 			  |
//-------------------------------------Cache Blocking Loop implementation End--------------------------------------

					long izEnd = MIN(bz+n3_Tblock, n3End);	
					long iyEnd = MIN(by+n2_Tblock, n2End);	
					long ixEnd = MIN(n1_Tblock, n1End-bx);
	
				if(bz==HALF_LENGTH){
					for(long i=0;i<bz; i++){
						for(long iy=by; iy<iyEnd; iy++) {
							//Compute the addresses for next x-loops
							tab = &tab_base[i*dimn1n2 + iy*n1 + bx];
	
							for(long ix=0;ix<ixEnd; ix++){
								tab[ix] = 0.f;
							}
						}
					}
				}	
				if(izEnd>=n3End){
					for(long i=n3End;i<n3; i++){
						for(long iy=by; iy<iyEnd; iy++) {
							//Compute the addresses for next x-loops
							tab = &tab_base[i*dimn1n2 + iy*n1 + bx];
	
							for(long ix=0;ix<ixEnd; ix++){
								tab[ix] = 0.f;
							}
						}
					}
				}
				if(by==HALF_LENGTH){
					for(long iz=bz; iz<izEnd; iz++) {
						for(long i=0;i<by; i++){
							//Compute the addresses for next x-loops
							tab = &tab_base[iz*dimn1n2 + i*n1 + bx];
	
							for(long ix=0;ix<ixEnd; ix++){
								tab[ix] = 0.f;
							}
						}
					}
				}
				if(iyEnd>=n2End){
					for(long iz=bz; iz<izEnd; iz++) {
						for(long i=n2End;i<n2; i++){
							//Compute the addresses for next x-loops
							tab = &tab_base[iz*dimn1n2 + i*n1 + bx];
	
							for(long ix=0;ix<ixEnd; ix++){
								tab[ix] = 0.f;
							}
						}
					}
				}

				if(bx==HALF_LENGTH){
					for(long iz=bz; iz<izEnd; iz++) {
						for(long iy=by; iy<iyEnd; iy++) {
							//Compute the addresses for next x-loops
							tab = &tab_base[iz*dimn1n2 + iy*n1 ];
	
							for(long ix=0;ix<HALF_LENGTH; ix++){
								tab[ix] = 0.f;
							}
						}
					}
				}
				if(ixEnd>=n1End){
					for(long iz=bz; iz<izEnd; iz++) {
						for(long iy=by; iy<iyEnd; iy++) {
							//Compute the addresses for next x-loops
							tab = &tab_base[iz*dimn1n2 + iy*n1 + ixEnd];
	
							for(long ix=0;ix<HALF_LENGTH; ix++){
								tab[ix] = 0.f;
							}
						}
					}
				}
					for(long iz=bz; iz<izEnd; iz++) {
						for(long iy=by; iy<iyEnd; iy++) {
							//Compute the addresses for next x-loops
							tab = &tab_base[iz*dimn1n2 + iy*n1 + bx];
	
							for(long ix=0;ix<ixEnd; ix++){
								tab[ix] = 0.f;
							}
						}
					}
				}
			}
		}
	
	}
}
#if defined(VERIFY_RESULTS)
#include <math.h>
void init_data(float *data, const long dimx, const long dimy, const long dimz)
{
  for(long iz=0; iz<dimz; iz++)
    for(long iy=0; iy<dimy; iy++)
      for(long ix=0; ix<dimx; ix++) {
	*data = (float)iz;
	++data;
      }
}


// naive and slow implementation
void reference_implementation(float *next, float *prev, float *coeff, 
		  float *vel,
		  const int n1, const int n2, const int n3, const int half_length){
  int n1n2 = n1*n2;
  
  for(int iz=0; iz<n3; iz++) {
    for(int iy=0; iy<n2; iy++) {
      for(int ix=0; ix<n1; ix++) {
	if( ix>=half_length && ix<(n1-half_length) && iy>=half_length && iy<(n2-half_length) && iz>=half_length && iz<(n3-half_length) ) {
	  float res = prev[iz*n1n2 + iy*n1 + ix]*coeff[0];
	  for(int ir=1; ir<=half_length; ir++) {
	    res += coeff[ir] * (prev[iz*n1n2 + iy*n1 + ix+ir] + prev[iz*n1n2 + iy*n1 + ix-ir]);	      // horizontal
	    res += coeff[ir] * (prev[iz*n1n2 + iy*n1 + ix+ir*n1] + prev[iz*n1n2 + iy*n1 + ix-ir*n1]);   // vertical
	    res += coeff[ir] * (prev[iz*n1n2 + iy*n1 + ix+ir*n1*n2] + prev[iz*n1n2 + iy*n1 + ix-ir*n1*n2]); // in front / behind
	  }
	  next[iz*n1n2 + iy*n1 +ix] = 2.0f* prev[iz*n1n2 + iy*n1 +ix] - next[iz*n1n2 + iy*n1 +ix] + res * vel[iz*n1n2 + iy*n1 +ix];
	}
      }
    }
  }
}

bool within_epsilon(float* output, float *reference, const long dimx, const long dimy, const long dimz, const long radius, const long zadjust=0, const float delta=0.0001f )
{
  bool retval = true;
  float abs_delta = fabsf(delta);
  for(long iz=0; iz<dimz; iz++) {
    for(long iy=0; iy<dimy; iy++) {
      for(long ix=0; ix<dimx; ix++) {
	if( ix>=radius && ix<(dimx-radius) && iy>=radius && iy<(dimy-radius) && iz>=radius && iz<(dimz-radius+zadjust) ) {
	  float difference = fabsf( *reference - *output);
	  if( difference > delta ) {
	    retval = false;
	    printf(" ERROR: (%d,%d,%d)\t%.2f instead of %.2f\n", ix,iy,iz, *output, *reference);
	    return false;
	  }
	}
	++output;
	++reference;
      }
    }
  }
  return retval;
}



#endif /* VERIFY_RESULTS */


#endif /*_TOOLS_INCLUDE */


