/*
  This file is provided under a BSD license.

  Copyright (c) 2005-2014 Intel Corporation. All rights reserved.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.
    * Neither the name of Intel Corporation nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#define TBB_PREVIEW_LOCAL_OBSERVER 1
#define TBB_PREVIEW_TASK_ARENA 1 /* needed until TBB version 4.3 */
#include <tbb/task_arena.h>
#include <tbb/task_scheduler_observer.h>
#include <errno.h>
#include <cstdlib>
#include <iostream>
#ifndef _WIN32
#include <sched.h>
#endif

namespace observers {

// Pin software threads to hardware resources. See the blog for more details:
// https://software.intel.com/en-us/blogs/2013/10/31/applying-intel-threading-building-blocks-observers-for-thread-affinity-on-intel
class pinning_observer: public tbb::task_scheduler_observer {
#ifdef _WIN32
typedef DWORD_PTR AFFINITY_MASK_TYPE;
#else
typedef cpu_set_t AFFINITY_MASK_TYPE;
#endif    
    bool pin_master;
    const int* cpu_ids;
    tbb::combinable<AFFINITY_MASK_TYPE> original_affinity;

public:
    pinning_observer( const int *_cpu_ids=NULL, bool _pin_master=true )
    : tbb::task_scheduler_observer(true), pin_master(_pin_master), cpu_ids(_cpu_ids) {}

    pinning_observer( tbb::task_arena &a, const int *_cpu_ids=NULL, bool _pin_master=false )
    : tbb::task_scheduler_observer(a), pin_master(_pin_master), cpu_ids(_cpu_ids) {}

    void set_map(const int *_cpu_ids) { cpu_ids = _cpu_ids; }

/*override*/ void on_scheduler_entry( bool is_worker ) {
        if ( !is_worker && !pin_master ) return;

        int thr_idx = tbb::task_arena::current_thread_index();
        int mapped_idx = cpu_ids ? cpu_ids[thr_idx] : thr_idx;

        AFFINITY_MASK_TYPE affinity, curr_affinity;

#ifndef _WIN32
        CPU_ZERO(&affinity); CPU_ZERO(&curr_affinity);
        CPU_SET(mapped_idx, &affinity);

        sched_getaffinity(0, sizeof(curr_affinity), &curr_affinity);
        int err = sched_setaffinity(0, sizeof(affinity), &affinity);
        if( err ) { perror("numa_node_to_cpus"); exit( EXIT_FAILURE ); }
#else
        affinity = ((AFFINITY_MASK_TYPE)1) << mapped_idx;
        curr_affinity = SetThreadAffinityMask(GetCurrentThread(), affinity);
#endif
        original_affinity.local() = curr_affinity;
    }

    /*override*/ void on_scheduler_exit( bool is_worker) {
        if ( !is_worker && !pin_master ) return;

        AFFINITY_MASK_TYPE affinity = original_affinity.local();
#ifndef _WIN32
        sched_setaffinity(0, sizeof(affinity), &affinity);
#else
        SetThreadAffinityMask(GetCurrentThread(), affinity);
#endif
    }
};
}
