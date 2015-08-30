// HyperThreadPhalanx.h
/***********************************************************************

 Copyright (c) 2013, James G. Dempsey

 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <assert.h>
#if !defined(__linux)
#include <intrin.h>
inline int __sync_fetch_and_add(int*Addend, int Value)
{
  return _InterlockedExchangeAdd((volatile int *)Addend, Value);
}
#endif
// types:
struct HyperThreadPhalanxThreadInfo_t
{
  int APICid;
  int PhysicalCore;
  int PhysicalHT;
  int LogicalCore;
  int LogicalHT;
};

struct HyperThreadPhalanx_t
{
  int isIntel;
  union {
  char ProcessorBrand[48];
  unsigned int ProcessorBrand_uint32[12];
  };
  int nHTsPerCore;      // hardware
  int nThreads;		// omp_get_num_threads() {in parallel region, no nesting}
  int nCores;		// number of core derived threfrom
  int nHTs;		// smallest number of HT's in mapped cores (logical HTs/core)
  struct HyperThreadPhalanxThreadInfo_t* ThreadInfo; // allocated to nThreads
};

// global variables:
extern struct HyperThreadPhalanx_t HyperThreadPhalanx;


// global thread private variables:
// myCore: 0-based logical core number (may be subset of physical cores and not necessarily including core(0))
// myHT:   0-based logical thread within core (may be subset of hw threads in core and not necessarily including hwThread(0))
#if defined(__linux)
extern __thread int myCore;
extern __thread int myHT;
#else
extern __declspec(thread) int myCore;
extern __declspec(thread) int myHT;
#endif

// functions:
int HyperThreadPhalanxInit();

