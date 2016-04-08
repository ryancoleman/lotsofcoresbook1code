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
 * ! Version 09
 * ! leonardo.borges@intel.com
 * ! cedric.andreolli@intel.com
 * !****************************************************************************/

#include "iso-3dfd.h"




void iso_3dfd_it(float *ptr_next_base,  float *ptr_prev_base,  float *ptr_vel_base,   float *coeff,
	      const TYPE_INTEGER n1, const TYPE_INTEGER n2, const TYPE_INTEGER n3, const TYPE_INTEGER num_threads,
	      const TYPE_INTEGER n1_Tblock, const TYPE_INTEGER n2_Tblock, const TYPE_INTEGER n3_Tblock){

	#pragma omp parallel OMP_N_THREADS
	{

		const TYPE_INTEGER dimn1n2 = n1*n2;
		const TYPE_INTEGER fullLength = 2*HALF_LENGTH;
	
		//One of the role of unrolling is to help the compiler to vectorize but also to help
		//it to fill the pipeline. In order to save some multiplications, we can save the 
		//values used to compute offsets.
		const TYPE_INTEGER vertical_1 = n1, vertical_2 = n1*2, vertical_3 = n1*3, vertical_4 = n1*4;
		const TYPE_INTEGER front_1 = dimn1n2, front_2 = dimn1n2*2, front_3 = dimn1n2*3, front_4 = dimn1n2*4;
		const float c0=coeff[0], c1=coeff[1], c2=coeff[2], c3=coeff[3], c4=coeff[4];
		//At this point, we must handle the stencil possible sizes.
	#if ( HALF_LENGTH == 8 )
		const TYPE_INTEGER vertical_5 = n1*5, vertical_6 = n1*6, vertical_7 = n1*7, vertical_8 = n1*8;
		const TYPE_INTEGER front_5 = dimn1n2*5, front_6 = dimn1n2*6, front_7 = dimn1n2*7, front_8 = dimn1n2*8;
		const float c5=coeff[5], c6=coeff[6], c7=coeff[7], c8=coeff[8];
	#endif

		const TYPE_INTEGER vertical[HALF_LENGTH] = {vertical_1, vertical_2, vertical_3, vertical_4
	#if ( HALF_LENGTH == 8 )
		, vertical_5, vertical_6, vertical_7, vertical_8
	#endif
		};
		const TYPE_INTEGER front[HALF_LENGTH] = {front_1, front_2, front_3, front_4
	#if ( HALF_LENGTH == 8 )
		, front_5, front_6, front_7, front_8
	#endif
		};

		TYPE_INTEGER n3End = n3 - HALF_LENGTH;
		TYPE_INTEGER n2End = n2 - HALF_LENGTH;
		TYPE_INTEGER n1End = n1 - HALF_LENGTH;


		register SIMD_TYPE yNextVec, tempVec;  

		register SIMD_TYPE beforeVec, afterVec;
		register SIMD_TYPE currentVec;
		SIMD_TYPE velVec, nextVec;                     
            	register SIMD_TYPE upVec;   

		register SIMD_TYPE divVec;

		register SIMD_TYPE coeffV;						
		float* ptr_next;
		register float* ptr_prev;
		float* ptr_vel;



		#if defined(AVX)
			SIMD_TYPE cVec, minusVec, plusVec, tempVec_1, yPrevVec, ySumVec, zNextVec, zPrevVec, zSumVec, tempSumVec, tempVec_3, tempVec_2;                     
		#endif

		register int offset;
		SIMD_TYPE twoVec = SET_VALUE_INTR(2.0f); 

		SIMD_TYPE coeffVec[HALF_LENGTH + 1];
		#pragma noprefetch
		for(TYPE_INTEGER i=0; i<=HALF_LENGTH; i++){
			GENERATE_COEFFICIENT_ARRAY_INTR(i)
		}
		
		
//-------------------------------------Cache Blocking Loop Implementation Begin----------------------------------
//												 		 |

	#ifdef BLOCK_X_Y_Z
		#pragma omp for OMP_SCHEDULE collapse(3)	
		for(TYPE_INTEGER bx=HALF_LENGTH; bx<n1End; bx+=n1_Tblock){	
			for(TYPE_INTEGER by=HALF_LENGTH; by<n2End; by+=n2_Tblock){	
				for(TYPE_INTEGER bz=HALF_LENGTH; bz<n3End; bz+=n3_Tblock){	
	#endif

	#ifdef BLOCK_X_Z_Y
		#pragma omp for OMP_SCHEDULE collapse(3)	
		for(TYPE_INTEGER bx=HALF_LENGTH; bx<n1End; bx+=n1_Tblock){	
			for(TYPE_INTEGER bz=HALF_LENGTH; bz<n3End; bz+=n3_Tblock){	
				for(TYPE_INTEGER by=HALF_LENGTH; by<n2End; by+=n2_Tblock){	
	#endif

	#ifdef BLOCK_Y_Z_X
	        #pragma omp for OMP_SCHEDULE collapse(1) 
	        for(TYPE_INTEGER by=HALF_LENGTH; by<n2End; by+=n2_Tblock){
	                for(TYPE_INTEGER bz=HALF_LENGTH; bz<n3End; bz+=n3_Tblock){       
	        		for(TYPE_INTEGER bx=HALF_LENGTH; bx<n1End; bx+=n1_Tblock){
	#endif

	#ifdef BLOCK_Y_X_Z
	        #pragma omp for OMP_SCHEDULE collapse(3) 
	        for(TYPE_INTEGER by=HALF_LENGTH; by<n2End; by+=n2_Tblock){
	                for(TYPE_INTEGER bx=HALF_LENGTH; bx<n1End; bx+=n1_Tblock){
	                	for(TYPE_INTEGER bz=HALF_LENGTH; bz<n3End; bz+=n3_Tblock){    
	#endif  

	#ifdef BLOCK_Z_X_Y
	        #pragma omp for OMP_SCHEDULE collapse(3) 
	        for(TYPE_INTEGER bz=HALF_LENGTH; bz<n3End; bz+=n3_Tblock){
	                for(TYPE_INTEGER bx=HALF_LENGTH; bx<n1End; bx+=n1_Tblock){
	        		for(TYPE_INTEGER by=HALF_LENGTH; by<n2End; by+=n2_Tblock){
	#endif 

	#ifdef BLOCK_Z_Y_X
	        #pragma omp for collapse(3) 
	       	for(TYPE_INTEGER bz=HALF_LENGTH; bz<n3End; bz+=n3_Tblock){
	        	for(TYPE_INTEGER by=HALF_LENGTH; by<n2End; by+=n2_Tblock){
	                	for(TYPE_INTEGER bx=HALF_LENGTH; bx<n1End; bx+=n1_Tblock){
	#endif 
//											 			  |
//-------------------------------------Cache Blocking Loop implementation End--------------------------------------

					TYPE_INTEGER izEnd = MIN(bz+n3_Tblock, n3End);	
					TYPE_INTEGER iyEnd = MIN(by+n2_Tblock, n2End);	
					TYPE_INTEGER ixEnd = MIN(n1_Tblock, n1End-bx);
					
					for(TYPE_INTEGER iz=bz; iz<izEnd; iz++) {
						for(TYPE_INTEGER iy=by; iy<iyEnd; iy++) {
							//Compute the addresses for next x-loops
							offset = iz*dimn1n2 + iy*n1 + bx;
							ptr_next = ptr_next_base + offset;
							ptr_prev = ptr_prev_base + offset;
							ptr_vel = ptr_vel_base + offset;

	
						//	#pragma vector nontemporal(ptr_prev)
							#pragma noprefetch
							#pragma ivdep
							for(TYPE_INTEGER ix=0;ix<ixEnd; ix+=SIMD_STEP){
							#if defined(MODEL_MIC)

								//----------------- x-dimension begin -----------------------
								SHIFT_MULT_INIT
									
								SHIFT_MULT_INTR(1, vertical_1, front_1, coeffVec[1])								
								SHIFT_MULT_INTR(2, vertical_2, front_2, coeffVec[2])									
								SHIFT_MULT_INTR(3, vertical_3, front_3, coeffVec[3])									
								SHIFT_MULT_INTR(4, vertical_4, front_4, coeffVec[4])		
							#if ( HALF_LENGTH == 8 )				
								SHIFT_MULT_INTR(5, vertical_5, front_5, coeffVec[5])									
								SHIFT_MULT_INTR(6, vertical_6, front_6, coeffVec[6])									
								SHIFT_MULT_INTR(7, vertical_7, front_7, coeffVec[7])									
								SHIFT_MULT_INTR(8, vertical_8, front_8, coeffVec[8])
							#endif
							#elif defined(AVX)
								SHIFT_MULT_INIT
									
								SHIFT_MULT_INTR(1)
								SHIFT_MULT_INTR(2)
								SHIFT_MULT_INTR(3)
								SHIFT_MULT_INTR(4)
							#if ( HALF_LENGTH == 8 )				
								SHIFT_MULT_INTR(5)
								SHIFT_MULT_INTR(6)
								SHIFT_MULT_INTR(7)
								SHIFT_MULT_INTR(8)
							#endif
									
								//----------------- x-dimension end --------------------------
								
								//----------------- y/z-dimension begin ----------------------
								MUL_COEFF_INTR(vertical_1, front_1, coeffVec[1])
								MUL_COEFF_INTR(vertical_2, front_2, coeffVec[2])
								MUL_COEFF_INTR(vertical_3, front_3, coeffVec[3])
								MUL_COEFF_INTR(vertical_4, front_4, coeffVec[4])
							#if ( HALF_LENGTH == 8 ) 
								MUL_COEFF_INTR(vertical_5, front_5, coeffVec[5])
								MUL_COEFF_INTR(vertical_6, front_6, coeffVec[6])
								MUL_COEFF_INTR(vertical_7, front_7, coeffVec[7])
								MUL_COEFF_INTR(vertical_8, front_8, coeffVec[8])
							#endif

								//----------------- y/z-dimension end ------------------------
							#endif //MODEL=AVX	
								//--------------  refreshing ptr_next begin ------------------
								REFRESH_NEXT_INTR
	

																//--------------  refreshing ptr_next end ---------------------
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


