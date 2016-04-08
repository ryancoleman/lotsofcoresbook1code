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


#ifndef _ISO_3DFD_INCLUDE
#define _ISO_3DFD_INCLUDE

#include <omp.h>
#include <stdio.h>
#include <immintrin.h>


/*** stencil compute behavior ***/


/**** Verify results after one ietration? *******/
//#define VERIFY_RESULTS

/*** verify if stencil half lenght is properly defined */
#if !( (HALF_LENGTH == 4) || (HALF_LENGTH == 8) ) 
#error "HALF_LENGTH must be defined (either 4 or 8)"
#endif

/**** Memory allocation ******/
//#define USE_MEMALIGN -- not in use
//#define LARGE_PAGES -- not in use
#define MASK_ALLOC_OFFSET(x) (x)
//#define MASK_ALLOC_OFFSET(x)  0
#define CACHELINE_BYTES   64


/* To offset halo and have arrays 64Byte aligned */
//#define ALIGN_HALO_FACTOR  (-HALF_LENGTH+1)
#define ALIGN_HALO_FACTOR  (-HALF_LENGTH )
/* #define  ALIGN_HALO_FACTOR  0     <-- no alignment */


/***** OpenMP schedule *********/
#define OMP_SCHEDULE schedule(dynamic)
//#define OMP_SCHEDULE schedule(guided)
//#define OMP_SCHEDULE schedule(static)



#define OMP_N_THREADS  num_threads(num_threads)

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define CEILING(X) (X-(int)(X) > 0 ? (int)(X+1) : (int)(X))

  #define P2 64 
#define P1 16

#define FINITE_ADD(ix, off) (ptr_prev[ix + off] + ptr_prev[ix - off])
#define C_PREFETCH1(addr, PREF_DIST) _mm_prefetch((char const*) ((addr) + (PREF_DIST)), _MM_HINT_T0);
#define C_PREFETCH2(addr, PREF_DIST) _mm_prefetch((char const*) ((addr) + (PREF_DIST)), _MM_HINT_T1);



#if defined(MODEL_MIC)
	#define TYPE_INTEGER int

	#define C_PREFETCH(ptr, v2, v1)                                                 \
                C_PREFETCH2(ptr, v2)                                                    \
                C_PREFETCH1(ptr, v1)

	#define SET_VALUE_INTR(v)							\
	 	_mm512_set1_ps(v)

	#define GENERATE_COEFFICIENT_ARRAY_INTR(i)                                      \
                C_PREFETCH(&coeff[i], 128, 16)                                           \
                coeffVec[i] = _mm512_set1_ps(coeff[i]);


        #define STORE_YZ_DIMENSION_IN_DIV_INTR                                          \
                C_PREFETCH2(div + ix, P2)                                               \
                C_PREFETCH1(div + ix, P1)                                               \
                _mm512_storenr_ps(div + ix, divVec);

	#define SHIFT_MULT_INTR(ind, vertical, front, coeff)							\
		coeffV = coeff;											\
               	yNextVec = (__m512)_mm512_alignr_epi32((__m512i)currentVec, (__m512i)beforeVec, SIMD_STEP-ind); \
 		__assume_aligned((void*)vertical,     CACHELINE_BYTES);                 \
                __assume_aligned((void*)&ptr_prev[ix+vertical],     CACHELINE_BYTES);   \
                C_PREFETCH(&ptr_prev[ix + vertical], 32, P1)                            \
      	    	divVec = _mm512_fmadd_ps(yNextVec, coeffV, divVec);                 				\
               	yNextVec = (__m512)_mm512_alignr_epi32((__m512i)afterVec, (__m512i)currentVec, ind); 		\
 		__assume_aligned((void*)&ptr_prev[ix-vertical],     CACHELINE_BYTES);   			\
                C_PREFETCH(&ptr_prev[ix - vertical], 32, P1)                            			\
		divVec = _mm512_fmadd_ps(yNextVec, coeffV, divVec);                 				\
                												\
                yNextVec = _mm512_load_ps(&ptr_prev[ix + vertical]);                    			\
       	    	divVec = _mm512_fmadd_ps(yNextVec, coeffV, divVec);                 				\
														\
                yNextVec = _mm512_load_ps(&ptr_prev[ix - vertical]);                    			\
 		__assume_aligned((void*)front,     CACHELINE_BYTES);                    			\
                __assume_aligned((void*)&ptr_prev[ix+front],     CACHELINE_BYTES);      			\
                C_PREFETCH(&ptr_prev[ix + front], 80, P1)                               			\
       	    	divVec = _mm512_fmadd_ps(yNextVec, coeffV, divVec);                 				\
                                                                                      				\
                yNextVec = _mm512_load_ps(&ptr_prev[ix + front]);                       			\
             	 __assume_aligned((void*)&ptr_prev[ix-front],     CACHELINE_BYTES);      			\
                C_PREFETCH(&ptr_prev[ix - front], 80, P1)                               			\
       	    	divVec = _mm512_fmadd_ps(yNextVec, coeffV, divVec);                 				\
														\
                yNextVec = _mm512_load_ps(&ptr_prev[ix - front]);                       			\
       	    	divVec = _mm512_fmadd_ps(yNextVec, coeffV, divVec);                 				\




	#define SHIFT_MULT_INIT								\
		__assume_aligned(&ptr_prev[ix-SIMD_STEP], CACHELINE_BYTES);		\
		C_PREFETCH(&ptr_prev[ix-SIMD_STEP], 256, 48)				\
		beforeVec = _mm512_load_ps(&ptr_prev[ix-SIMD_STEP]);			\
		currentVec = _mm512_load_ps(&ptr_prev[ix]);				\
		afterVec = _mm512_load_ps(&ptr_prev[ix+SIMD_STEP]);			\
		divVec = _mm512_mul_ps(currentVec, coeffVec[0]);			\
											
	
	


	#define REFRESH_NEXT_INTR									\
		C_PREFETCH2(&ptr_next[ix], 64)								\
		C_PREFETCH1(&ptr_next[ix], 16)								\
		nextVec = _mm512_load_ps(&ptr_next[ix]); 						\
		C_PREFETCH2(&ptr_vel[ix], 64)								\
		C_PREFETCH1(&ptr_vel[ix], 16)								\
		velVec = _mm512_load_ps(&ptr_vel[ix]);							\
		nextVec = _mm512_fmadd_ps(divVec, velVec, _mm512_fmsub_ps(currentVec, twoVec, nextVec));\
		_mm512_storenrngo_ps(&ptr_next[ix], nextVec);



#elif defined(AVX)
	#define TYPE_INTEGER long

        #define C_PREFETCH(ptr, v2, v1)                                         


        #define SET_ZERO_INTR _mm256_setzero_ps()
        

	#define SET_VALUE_INTR(v) _mm256_set1_ps(v)


        #define GENERATE_COEFFICIENT_ARRAY_INTR(i)                                      \
                coeffVec[i] = _mm256_set1_ps(coeff[i]);


        #define STORE_YZ_DIMENSION_IN_DIV_INTR                                          \
                _mm256_store_ps(div + ix, divVec);


        #define MUL_COEFF_INTR(vertical, front, coeff)                                  \
                __assume_aligned((void*)vertical,     CACHELINE_BYTES);                 \
                __assume_aligned((void*)&ptr_prev[ix+vertical],     CACHELINE_BYTES);   \
                yNextVec = _mm256_load_ps(&ptr_prev[ix + vertical]);                    \
                                                                                        \
                __assume_aligned((void*)&ptr_prev[ix-vertical],     CACHELINE_BYTES);   \
                yPrevVec = _mm256_load_ps(&ptr_prev[ix - vertical]);                    \
                ySumVec = _mm256_add_ps(yNextVec, yPrevVec);                            \
                                                                                        \
                __assume_aligned((void*)front,     CACHELINE_BYTES);                    \
                __assume_aligned((void*)&ptr_prev[ix+front],     CACHELINE_BYTES);      \
                zNextVec = _mm256_load_ps(&ptr_prev[ix + front]);                       \
                                                                                        \
                __assume_aligned((void*)&ptr_prev[ix-front],     CACHELINE_BYTES);      \
                zPrevVec = _mm256_load_ps(&ptr_prev[ix - front]);                       \
                                                                                        \
                zSumVec = _mm256_add_ps(zNextVec, zPrevVec);                            \
                                                                                        \
                tempSumVec = _mm256_add_ps(ySumVec, zSumVec);                           \
                                                                                        \
                zSumVec = _mm256_mul_ps(tempSumVec, coeff);                             \
                divVec = _mm256_add_ps(zSumVec, divVec);


	
	#define SHIFT_MULT_INIT								\
		__assume_aligned(&ptr_prev[ix-SIMD_STEP], CACHELINE_BYTES);		\
											\
		beforeVec = _mm256_load_ps(&ptr_prev[ix-SIMD_STEP]);			\
		currentVec = _mm256_load_ps(&ptr_prev[ix]);				\
		afterVec = _mm256_load_ps(&ptr_prev[ix+SIMD_STEP]);			\
		divVec = _mm256_mul_ps(currentVec, coeffVec[0]);			
											
	
	
	#define SHIFT_MULT_INTR(ind)					\
		cVec = coeffVec[ind];  					\
				          				\
               	minusVec = _mm256_load_ps(&ptr_prev[ix - ind]); 	\
               	plusVec = _mm256_load_ps(&ptr_prev[ix +ind]); 		\
       	    	upVec = _mm256_add_ps(minusVec, plusVec);               \
		tempVec_1 = _mm256_mul_ps(upVec, cVec);			\
		divVec = _mm256_add_ps(tempVec_1, divVec);    		\             					




	#define REFRESH_NEXT_INTR									\
		velVec = _mm256_load_ps(&ptr_vel[ix]);							\
		tempVec_3 = _mm256_mul_ps(divVec, velVec);/*div*vel -> temp_3*/				\
                        		 								\
		__assume_aligned((void*)&ptr_next[ix],     CACHELINE_BYTES);				\
		nextVec = _mm256_load_ps(&ptr_next[ix]); 						\
		tempVec_1 = _mm256_mul_ps(twoVec, currentVec);/*2*prev -> temp_1*/			\
		tempVec_2 = _mm256_sub_ps(tempVec_1, nextVec);/* 2*prev - next -> temp_2*/		\
		nextVec = _mm256_add_ps(tempVec_3, tempVec_2);/*2*prev - next + div * vel -> next*/ 	\
		_mm256_store_ps(&ptr_next[ix], nextVec);

#endif 

void iso_3dfd_stencil_BLK(float *ptr_next,  float *ptr_prev,  float *ptr_vel, float *coeff,
			  const TYPE_INTEGER i1_idx, const TYPE_INTEGER i2_idx, const TYPE_INTEGER i3_idx,
			  const TYPE_INTEGER n1,     const TYPE_INTEGER n2,     const TYPE_INTEGER n3,
			  const TYPE_INTEGER n1_Tb,  const TYPE_INTEGER n2_Tb,  const TYPE_INTEGER n3_Tb);

void iso_3dfd(float *next,  float *prev,  float *vel,   float *coeff,
	      const TYPE_INTEGER n1, const TYPE_INTEGER n2, const TYPE_INTEGER n3, const TYPE_INTEGER num_threads, const TYPE_INTEGER nreps,
	      const TYPE_INTEGER n1_Tblock, const TYPE_INTEGER n2_Tblock, const TYPE_INTEGER n3_Tblock);

#define PRAGMA_LONG_LOOP _Pragma("loop_count(100000)")

#endif /* _ISO_3DFD_INCLUDE */
