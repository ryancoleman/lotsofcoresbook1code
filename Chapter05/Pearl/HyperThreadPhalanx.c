// HyperThreadPhalanx.c
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

#include "HyperThreadPhalanx.h"

struct HyperThreadPhalanx_t HyperThreadPhalanx;

#if defined(__linux)
__thread int myCore = -1;	// logical core (may be subset of physical cores and not necessarily core(0))
__thread int myHT = -1;	// logical thread (may be subset of hw threads in core and not necessarily hwThread(0) in core)
#else
__declspec(thread) int myCore = -1;	// logical core (may be subset of physical cores and not necessarily core(0))
__declspec(thread) int myHT = -1;	// logical thread (may be subset of hw threads in core and not necessarily hwThread(0) in core)
#endif

#if !defined(BOOL)
#define BOOL int
#endif

#if !defined(TRUE)
#define TRUE 1
#endif
#if !defined(FALSE)
#define FALSE 0
#endif

#if defined(__MIC__)
#define WAIT_A_BIT _mm_delay_32(10)
#else
#define WAIT_A_BIT _mm_pause();
#endif

void __cpuidEX(int cpuinfo[4], int func_a, int func_c)
{
	int eax, ebx, ecx, edx;
	__asm__ __volatile__ ("cpuid":\
	"=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx) : "a" (func_a), "c" (func_c));
	cpuinfo[0] = eax;
	cpuinfo[1] = ebx;
	cpuinfo[2] = ecx;
	cpuinfo[3] = edx;
} // void __cpuidEX(int cpuinfo[4], int func_a, int func_c)

void InitProcessor()
{
  unsigned int CPUinfo[4];
  __cpuid(CPUinfo, 0);	// This code requires at least support of CPUID
  HyperThreadPhalanx.ProcessorBrand_uint32[0] = CPUinfo[1];
  HyperThreadPhalanx.ProcessorBrand_uint32[1] = CPUinfo[3];	// note order different
  HyperThreadPhalanx.ProcessorBrand_uint32[2] = CPUinfo[2];
  HyperThreadPhalanx.ProcessorBrand_uint32[3] = 0;
  HyperThreadPhalanx.isIntel = (strcmp(HyperThreadPhalanx.ProcessorBrand, "GenuineIntel") == 0);
}

int HyperThreadPhalanxInit()
{
  InitProcessor();
  if(!HyperThreadPhalanx.isIntel)
  {
    printf("Not Intel processor. Add code to handle this processor.\n");
    return 1;
  }
  if(omp_in_parallel())
  {
    printf("HyperThreadPhalanxInit() must be called from outside parallel .\n");
    return 2;
  }

  BOOL GlobalOverSubscription = FALSE;
  BOOL GlobalUnderSubscription = FALSE;
  volatile int OverSubscriptionThreadCount;
#pragma omp parallel
  {
    int nThreads = omp_get_num_threads();	// use omp_get_num_threads() NOT omp_get_max_threads()
    int iThread = omp_get_thread_num();
    unsigned int CPUinfo[4];

#pragma omp master
    {
      HyperThreadPhalanx.nThreads = nThreads;
      HyperThreadPhalanx.ThreadInfo = malloc(nThreads * sizeof(struct HyperThreadPhalanxThreadInfo_t));
      __cpuidEX(CPUinfo, 4, 0);
      HyperThreadPhalanx.nHTsPerCore = ((CPUinfo[0] >> 14) & 0x3F) + 1;
      HyperThreadPhalanx.nHTs = HyperThreadPhalanx.nHTsPerCore; // set logical HT's per core to physical
    }
#pragma omp barrier
    // OpenMP implimentation detail:
    //
    // It appears that on some systems (Host in particular) that KMP_AFFINITY=compact or scatter
    // binds a software thread to a core's first level cache, and not to a logical processor.
    // IOW the software thread is free to migrate amongst the hardware threads within the core.
    // Whereas on other systems the software thread is bound to a single hardware thread.
    // To overcome this we will put the program into a short-lived compute session that obtains
    // one of the mappings of software thread to hardware thread. The HyperThread Phalanx will
    // not be sensitive to logically mapped HTs migrating within the original core placement.
    // When the system is well behaved the iMapAttempt loop will iterate once.
    BOOL PrivateOverSubscription = FALSE;
#pragma omp barrier
#pragma omp master
      OverSubscriptionThreadCount = 0;
#pragma omp barrier
      __sync_fetch_and_add((int*)&OverSubscriptionThreadCount, 1);
      // sit in minor spinwait until all nThreads have incrimented the count
      // if the user did not over-subscribe, 
      // all software threads will be running on different hardware threads
      while(OverSubscriptionThreadCount < nThreads)
        WAIT_A_BIT;
      // master region finished, see if allocation succeded
      if(HyperThreadPhalanx.ThreadInfo)
      {
        __cpuidEX(CPUinfo, 1, 0); // get features
        if(CPUinfo[2] & (1 << 21))
        {
          // processor has x2APIC
          __cpuidEX(CPUinfo, 0x0B, 0);
	  // get thread's APICid
          HyperThreadPhalanx.ThreadInfo[iThread].APICid = CPUinfo[3];
        }
        else
        {
          // older processor without x2APIC
          __cpuidEX(CPUinfo, 1, 0);
	  // get thread's APICid
          HyperThreadPhalanx.ThreadInfo[iThread].APICid = (CPUinfo[1] >> 24) & 0xFF;
        }
        // Use thread's APICid to determine physical core and physical HT number within core
        HyperThreadPhalanx.ThreadInfo[iThread].PhysicalCore = HyperThreadPhalanx.ThreadInfo[iThread].APICid / HyperThreadPhalanx.nHTsPerCore;
        HyperThreadPhalanx.ThreadInfo[iThread].PhysicalHT = HyperThreadPhalanx.ThreadInfo[iThread].APICid % HyperThreadPhalanx.nHTsPerCore;
        HyperThreadPhalanx.ThreadInfo[iThread].LogicalCore = -1;	// for now indicate not assigned
        HyperThreadPhalanx.ThreadInfo[iThread].LogicalHT = -1;	// for now indicate not assigned
      }
#pragma omp barrier
      // At this point, all the HyperThreadPhalanx.ThreadInfo[iThread].APICid, PhysicalCore and PhysicalHT have been filled-in
      // However, the logical core number may differ from physical core number
      // no different than OpenMP thread number differing from logical processor number
      // The logical core numbers are 0-based, without gaps
#pragma omp master
      {
        int NextLogicalCore = 0;
        for(;;)
        {
          int iLowest = -1; // none found
          for(int i = 0; i < HyperThreadPhalanx.nThreads; ++i)
          {
            // see if unassigned core
            if(HyperThreadPhalanx.ThreadInfo[i].LogicalCore == -1)
            {
              if(iLowest < 0)
              {
                // first unassigned is lowest
                iLowest = i;
              }
              else
              {
                if(HyperThreadPhalanx.ThreadInfo[i].APICid < HyperThreadPhalanx.ThreadInfo[iLowest].APICid)
              	  iLowest = i;	// new lowest
              }
            } // if(HyperThreadPhalanx.ThreadInfo[i].LogicalCore < 0)
          } // for(int i = 0; i < HyperThreadPhalanx.nThreads; ++i)
          if(iLowest < 0)
            break;
          // able to use core
          int NextLogicalHT = 0;
          for(int i = 0; i < HyperThreadPhalanx.nThreads; ++i)
          {
            if(HyperThreadPhalanx.ThreadInfo[i].PhysicalCore == HyperThreadPhalanx.ThreadInfo[iLowest].PhysicalCore)
            {
              HyperThreadPhalanx.ThreadInfo[i].LogicalCore = NextLogicalCore; // mark as unavailable with NextLogicalCore
              HyperThreadPhalanx.ThreadInfo[i].LogicalHT = NextLogicalHT++;   // assign NextLogicalHT and advance
              if(NextLogicalHT > HyperThreadPhalanx.nHTs)
                GlobalOverSubscription = TRUE;	// only when observed, not always
            }
          } // for(int i = 0; i < HyperThreadPhalanx.nThreads; ++i)
          ++NextLogicalCore;
          if(NextLogicalHT < HyperThreadPhalanx.nHTs)
            GlobalUnderSubscription = TRUE;
        } // for(;;)
        HyperThreadPhalanx.nCores = NextLogicalCore;
      } // omp master
#pragma omp barrier
      // master is finished
      myCore = HyperThreadPhalanx.ThreadInfo[iThread].LogicalCore;
      myHT = HyperThreadPhalanx.ThreadInfo[iThread].LogicalHT;
  } // omp parallel
  if(GlobalOverSubscription || GlobalUnderSubscription)
  {
    printf("Oversubscription and/or Undersubscription of threads\n");
    printf("SW threads assigned to HW thread\n"); 
    for(int iThread = 0; iThread < HyperThreadPhalanx.nThreads; ++iThread)
    {
      printf("iThread = %d, LogicalCore = %d, LogicalHT = %d\n",
		iThread,
                HyperThreadPhalanx.ThreadInfo[iThread].LogicalCore,
                HyperThreadPhalanx.ThreadInfo[iThread].LogicalHT);
    }
    return 4;
  }
  return 0;
} // void HyperThreadPhalanxInit()

