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

#include "iso-3dfd.h"

#include <omp.h>
#include <stdio.h>


void iso_3dfd_it(float *ptr_next_base,  float *ptr_prev_base,  float *ptr_vel_base, float *coeff,
	      const TYPE_INTEGER n1, const TYPE_INTEGER n2, const TYPE_INTEGER n3, const TYPE_INTEGER num_threads,
	      const TYPE_INTEGER n1_Tblock, const TYPE_INTEGER n2_Tblock, const TYPE_INTEGER n3_Tblock){
		#pragma omp parallel OMP_N_THREADS default(shared)
	{
		TYPE_INTEGER dimn1n2 = n1*n2;
		TYPE_INTEGER n3End = n3 - HALF_LENGTH;
		TYPE_INTEGER n2End = n2 - HALF_LENGTH;
		TYPE_INTEGER n1End = n1 - HALF_LENGTH;

		float* ptr_next;
		float* ptr_prev;
		float* ptr_vel;
		float value;
		TYPE_INTEGER izEnd;	
		TYPE_INTEGER iyEnd;	
		TYPE_INTEGER ixEnd;
		
		const TYPE_INTEGER vertical_1 = n1, vertical_2 = n1*2, vertical_3 = n1*3, vertical_4 = n1*4;
        	const TYPE_INTEGER front_1 = dimn1n2, front_2 = dimn1n2*2, front_3 = dimn1n2*3, front_4 = dimn1n2*4;
        	const float c0=coeff[0], c1=coeff[1], c2=coeff[2], c3=coeff[3], c4=coeff[4];
        	//At this poTYPE_INTEGER, we must handle the stencil possible sizes.
        #if ( HALF_LENGTH == 8 )
		const TYPE_INTEGER vertical_5 = n1*5, vertical_6 = n1*6, vertical_7 = n1*7, vertical_8 = n1*8;
		const TYPE_INTEGER front_5 = dimn1n2*5, front_6 = dimn1n2*6, front_7 = dimn1n2*7, front_8 = dimn1n2*8;
		const float c5=coeff[5], c6=coeff[6], c7=coeff[7], c8=coeff[8];
	#endif
		__assume_aligned((void*)vertical_1, CACHELINE_BYTES);
		__assume_aligned((void*)vertical_2, CACHELINE_BYTES);
		__assume_aligned((void*)vertical_3, CACHELINE_BYTES);
		__assume_aligned((void*)vertical_4, CACHELINE_BYTES);
		__assume_aligned((void*)front_1, CACHELINE_BYTES);
		__assume_aligned((void*)front_2, CACHELINE_BYTES);
		__assume_aligned((void*)front_3, CACHELINE_BYTES);
		__assume_aligned((void*)front_4, CACHELINE_BYTES);
		//Handle all size of stencil
	#if( HALF_LENGTH == 8 )
		__assume_aligned((void*)vertical_5, CACHELINE_BYTES);
		__assume_aligned((void*)vertical_6, CACHELINE_BYTES);
		__assume_aligned((void*)vertical_7, CACHELINE_BYTES);
		__assume_aligned((void*)vertical_8, CACHELINE_BYTES);
		__assume_aligned((void*)front_5, CACHELINE_BYTES);
		__assume_aligned((void*)front_6, CACHELINE_BYTES);
		__assume_aligned((void*)front_7, CACHELINE_BYTES);
		__assume_aligned((void*)front_8, CACHELINE_BYTES);
	#endif

		__declspec(align(CACHELINE_BYTES)) float div[n1_Tblock];

		#pragma omp for OMP_SCHEDULE collapse(3)	
		for(TYPE_INTEGER bz=HALF_LENGTH; bz<n3End; bz+=n3_Tblock){	
			for(TYPE_INTEGER by=HALF_LENGTH; by<n2End; by+=n2_Tblock){	
				for(TYPE_INTEGER bx=HALF_LENGTH; bx<n1End; bx+=n1_Tblock){	
					izEnd = MIN(bz+n3_Tblock, n3End);	
					iyEnd = MIN(by+n2_Tblock, n2End);	
					ixEnd = MIN(n1_Tblock, n1End-bx);
	
					for(TYPE_INTEGER iz=bz; iz<izEnd; iz++) {
						for(TYPE_INTEGER iy=by; iy<iyEnd; iy++) {
        				        	ptr_next = &ptr_next_base[iz*dimn1n2 + iy*n1 + bx];
							ptr_prev = &ptr_prev_base[iz*dimn1n2 + iy*n1 + bx];
							ptr_vel = &ptr_vel_base[iz*dimn1n2 + iy*n1 + bx];
						
							__assume_aligned(ptr_next, CACHELINE_BYTES);
							__assume_aligned(ptr_prev, CACHELINE_BYTES);
							__assume_aligned(ptr_vel, CACHELINE_BYTES);
							
							#pragma prefetch ptr_prev:1:2 
							#pragma prefetch ptr_prev:0:1
							#pragma vector nontemporal(value)
							#pragma ivdep
							for(TYPE_INTEGER ix=0; ix<ixEnd; ix++) {
								value = ptr_prev[ix]*c0
									+ c1 *   (FINITE_ADD(ix, 1)
										+ FINITE_ADD(ix, vertical_1) 
										+ FINITE_ADD(ix, front_1))
									+ c2 * 	 (FINITE_ADD(ix, 2)
										+ FINITE_ADD(ix, vertical_2)
										+ FINITE_ADD(ix, front_2))
									+ c3 * 	 (FINITE_ADD(ix, 3)
										+ FINITE_ADD(ix, vertical_3)
										+ FINITE_ADD(ix, front_3))
									+ c4 * 	 (FINITE_ADD(ix, 4)
										+ FINITE_ADD(ix, vertical_4)
										+ FINITE_ADD(ix, front_4))
							#if( HALF_LENGTH == 8)
									+ c5 * 	 (FINITE_ADD(ix, 5)
										+ FINITE_ADD(ix, vertical_5)
										+ FINITE_ADD(ix, front_5))
									+ c6 * 	 (FINITE_ADD(ix, 6)
										+ FINITE_ADD(ix, vertical_6)
										+ FINITE_ADD(ix, front_6))
									+ c7 * 	 (FINITE_ADD(ix, 7)
										+ FINITE_ADD(ix, vertical_7)
										+ FINITE_ADD(ix, front_7))
									+ c8 * 	 (FINITE_ADD(ix, 8)
										+ FINITE_ADD(ix, vertical_8)
										+ FINITE_ADD(ix, front_8))
							#endif
								;
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
	      const TYPE_INTEGER n1, const TYPE_INTEGER n2, const TYPE_INTEGER n3, const TYPE_INTEGER num_threads, const TYPE_INTEGER nreps,
	      const TYPE_INTEGER n1_Tblock, const TYPE_INTEGER n2_Tblock, const TYPE_INTEGER n3_Tblock){
 
  for(TYPE_INTEGER it=0; it<nreps; it+=2){

    iso_3dfd_it(ptr_next, ptr_prev, ptr_vel, coeff,
		n1, n2, n3, num_threads, n1_Tblock, n2_Tblock, n3_Tblock);

    // here's where boundary conditions+halo exchanges happen

    // Swap previous & next between iterations
    iso_3dfd_it(ptr_prev, ptr_next, ptr_vel, coeff,
		n1, n2, n3, num_threads, n1_Tblock, n2_Tblock, n3_Tblock);


  } // time loop
}


