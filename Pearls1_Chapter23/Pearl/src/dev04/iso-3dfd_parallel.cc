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
 * ! Version 04
 * ! leonardo.borges@intel.com
 * ! cedric.andreolli@intel.com
 * !****************************************************************************/

#include "iso-3dfd.h"

#include <omp.h>
#include <stdio.h>


void iso_3dfd_it(float *ptr_next_base,  float *ptr_prev_base,  float *ptr_vel_base, float *coeff,
	      const int n1, const int n2, const int n3, const int num_threads,
	      const int n1_Tblock, const int n2_Tblock, const int n3_Tblock){
	int dimn1n2 = n1*n2;
	int n3End = n3 - HALF_LENGTH;
	int n2End = n2 - HALF_LENGTH;
	int n1End = n1 - HALF_LENGTH;

	#pragma omp parallel OMP_N_THREADS default(shared)
	{
		float* ptr_next;
		float* ptr_prev;
		float* ptr_vel;
		float value;
		int izEnd;	
		int iyEnd;	
		int ixEnd;
	
		#pragma omp for OMP_SCHEDULE collapse(3)	
		for(int bz=HALF_LENGTH; bz<n3End; bz+=n3_Tblock){	
			for(int by=HALF_LENGTH; by<n2End; by+=n2_Tblock){	
				for(int bx=HALF_LENGTH; bx<n1End; bx+=n1_Tblock){	
					izEnd = MIN(bz+n3_Tblock, n3End);	
					iyEnd = MIN(by+n2_Tblock, n2End);	
					ixEnd = MIN(n1_Tblock, n1End-bx);
	
					for(int iz=bz; iz<izEnd; iz++) {
						for(int iy=by; iy<iyEnd; iy++) {
        				        	ptr_next = &ptr_next_base[iz*dimn1n2 + iy*n1 + bx];
							ptr_prev = &ptr_prev_base[iz*dimn1n2 + iy*n1 + bx];
							ptr_vel = &ptr_vel_base[iz*dimn1n2 + iy*n1 + bx];
							
							#pragma ivdep
							for(int ix=0; ix<ixEnd; ix++) {
								value = 0.0;
								value += ptr_prev[ix]*coeff[0];
							//Unrolling on CPU allows to vectorize
							#if defined(MODEL_CPU)
								#pragma unroll(16)
							#endif
								#pragma ivdep
								for(int ir=1; ir<=HALF_LENGTH; ir++) {
									value += coeff[ir] * (ptr_prev[ix + ir] + ptr_prev[ix - ir]);// horizontal
									value += coeff[ir] * (ptr_prev[ix + ir*n1] + ptr_prev[ix - ir*n1]);// vertical
									value += coeff[ir] * (ptr_prev[ix + ir*dimn1n2] + ptr_prev[ix - ir*dimn1n2]); // in front / behind
								}
								ptr_next[ix] = 2.0f* ptr_prev[ix] - ptr_next[ix] + value*ptr_vel[ix];
							}
						}
					}
				}
			}
		}
	}
	
}


/***************************************************************
 *
 * iso_3dfd: Creates a 3D thread blocking using
 *                 n1_Tblock x n2_Tblock x n3_Tblock
 *               blocks
 *
 ***************************************************************/

void iso_3dfd(float *ptr_next,  float *ptr_prev,  float *ptr_vel,   float *coeff,
	      const int n1, const int n2, const int n3, const int num_threads, const int nreps,
	      const int n1_Tblock, const int n2_Tblock, const int n3_Tblock){
 
  for(int it=0; it<nreps; it+=2){

    iso_3dfd_it(ptr_next, ptr_prev, ptr_vel, coeff,
		n1, n2, n3, num_threads, n1_Tblock, n2_Tblock, n3_Tblock);

    // here's where boundary conditions+halo exchanges happen

    // Swap previous & next between iterations
    iso_3dfd_it(ptr_prev, ptr_next, ptr_vel, coeff,
		n1, n2, n3, num_threads, n1_Tblock, n2_Tblock, n3_Tblock);


  } // time loop
}


