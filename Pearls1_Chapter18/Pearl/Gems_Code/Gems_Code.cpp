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

#include <random>
#include <chrono>
#include <memory>
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>

#include <tbb/compat/thread>
#include <tbb/combinable.h>
#include <tbb/atomic.h>
#include <tbb/tbb_thread.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/pipeline.h>
#include <tbb/cache_aligned_allocator.h>

#include "arenas.h"

#include "../ittnotify/ittnotify.h"

#ifdef __CILK_BASELINE__
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <atomic>
#endif

#ifdef __OMP_BASELINE__
#include <omp.h>
#endif

#define DATA_SIZE  (32*1024*1024)
#define BLOCK_SIZE (4*1024)
#define ITERATIONS 10

#ifdef __MIC__
    #ifndef MIC_TEAMS
        #define MIC_TEAMS 20
    #endif
    #ifndef _NUM_TASKS
        #define _NUM_TASKS      20
    #endif
#else
    #ifndef _NUM_TASKS
        #define _NUM_TASKS      60
    #endif
#endif
typedef std::vector<double, tbb::cache_aligned_allocator<double> > data_vector_t;
typedef std::vector<std::chrono::duration<double> >                duration_array_t;

static __itt_domain* ittDomain = __itt_domain_create("ArenaSample.Domain");
static __itt_string_handle* ittMemory = __itt_string_handle_create("MemoryHandling");
static __itt_string_handle* ittInput = __itt_string_handle_create("InputGeneration");
static __itt_string_handle* ittLoop = __itt_string_handle_create("Loop");
static __itt_string_handle* ittCoreFunction = __itt_string_handle_create("CoreFunction");

// Create all worker threads
void initialize_threading()
{
#if defined(__OMP_BASELINE__)
    // Create parallel region -> initilized all worker threads
    #pragma omp parallel
    {}
#elif defined(__CILK_BASELINE__)
    std::atomic<int> barrier;
    barrier = 0;
    const int num_workers = __cilkrts_get_nworkers();
    cilk_for(int i=0; i<num_workers; i++)
    {
        barrier++;
        while(barrier < num_workers);
    }
#elif defined(__TBB_BASELINE__)
    tbb::atomic<int> barrier;
    barrier = 0;
    const int num_workers = tbb::task_scheduler_init::default_num_threads();
    tbb::parallel_for(0,num_workers,[&](int i)
    {
        barrier++;
        while(barrier < num_workers);
    });
#else
    std::cerr << "Invalid function call" << std::endl;
    exit(-1);
#endif
}

// Measeure execution time of specified functor block
template <typename F> std::chrono::duration<double> time(const F& _f) {
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    _f();
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();

    return end - start;
}

void fill_input(double* dst, int begin, int end)
{
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0., 1.);
    auto rand_val = std::bind ( distribution, generator );
    
    for( int j=begin;j<end;++j) {
        dst[j] = rand_val();
    }
}

void compute_sqrt_block(const double* src, double* dst, int begin, int end, int size)
{
    __itt_task_begin(ittDomain, __itt_null, __itt_null, ittCoreFunction);

#if (__INTEL_COMPILER>=1500) || ((__GNUC__*10000+__GNUC_MINOR__*100+__GNUC_PATCHLEVEL__)>=40900)
    #pragma omp simd aligned(src, dst: 64)
#elif defined(__INTEL_COMPILER)
    #pragma simd aligned
#endif
    for(int i=begin; i<end; ++i)
    {
        dst[i] = sqrt(src[i]);
    }
    __itt_task_end(ittDomain);
}

void compute_task() {
    data_vector_t src_vect(DATA_SIZE);
    data_vector_t dst_vect(DATA_SIZE);

    double* src = &src_vect[0];
    double* dst = &dst_vect[0];

    // Fill source buffer with random data
    // Insure that memory is paged in for output buffer as weel
    __itt_task_begin(ittDomain, __itt_null, __itt_null, ittInput);

#if defined(__OMP_BASELINE__)
    #pragma omp parallel for
    for(int i=0; i<DATA_SIZE; i+=BLOCK_SIZE) {
        int begin = i, end = std::min(i+BLOCK_SIZE, DATA_SIZE);
        #pragma noinline
        fill_input(src, begin, end);
    }
    __itt_task_end(ittDomain);

    __itt_task_begin(ittDomain, __itt_null, __itt_null, ittLoop);
    for(int t=0; t<ITERATIONS; ++t) {
        #pragma omp parallel for
        for(int i=0; i<DATA_SIZE; i+=BLOCK_SIZE) {
            int begin = i, end = std::min(i+BLOCK_SIZE, DATA_SIZE);
            #pragma noinline
            compute_sqrt_block(src, dst, begin, end, DATA_SIZE);
        }
    }
#elif defined(__CILK_BASELINE__)
    cilk_for(int i=0; i<DATA_SIZE; i+=BLOCK_SIZE) {
        int begin = i, end = std::min(i+BLOCK_SIZE, DATA_SIZE);
        #pragma noinline
        fill_input(src, begin, end);
    }
    __itt_task_end(ittDomain);

    __itt_task_begin(ittDomain, __itt_null, __itt_null, ittLoop);
    for(int t=0; t<ITERATIONS; ++t) {
        cilk_for(int i=0; i<DATA_SIZE; i+=BLOCK_SIZE) {
            int begin = i, end = std::min(i+BLOCK_SIZE, DATA_SIZE);
            #pragma noinline
            compute_sqrt_block(src, dst, begin, end, DATA_SIZE);
        }
    }
#else
    tbb::parallel_for(0, DATA_SIZE, BLOCK_SIZE,
        [&](int i) {
            int begin = i, end = std::min(i+BLOCK_SIZE, DATA_SIZE);
            #pragma noinline
            fill_input(src, begin, end);
        });
    __itt_task_end(ittDomain);

    tbb::affinity_partitioner partitioner;
    __itt_task_begin(ittDomain, __itt_null, __itt_null, ittLoop);
    for(int t=0; t<ITERATIONS; ++t) {
        tbb::parallel_for(0, DATA_SIZE, BLOCK_SIZE,
            [&](int i)
            {
                int begin = i, end = std::min(i+BLOCK_SIZE, DATA_SIZE);
                #pragma noinline
                compute_sqrt_block(src, dst, begin, end, DATA_SIZE);
            }, partitioner);
    }
#endif
    __itt_task_end(ittDomain);
}

void execute_baseline(duration_array_t& task_duration) {
#if defined(__OMP_BASELINE__)
    std::cerr << "Executing OpenMP baseline" << std::endl;
    omp_set_dynamic(1);
    omp_set_nested(1);
    setenv("KMP_AFFINITY","granularity=fine,scatter",1);
    #pragma omp parallel
    {
        #pragma omp master
        {
            // though this is not optimal to use for-loop+task, it emulates adding tasks from different threads or/and modules
            for(int task=0; task<_NUM_TASKS; ++task) {
                #pragma omp task firstprivate(task)
                {
                    task_duration[task] = time([&] { // Start time measurment
                        compute_task();
                    }); // end of measure
                }
            }
        }
    }
#elif defined(__CILK_BASELINE__)
    std::cerr << "Executing Cilk baseline" << std::endl;
    // though this is not optimal to use for-loop+cilk_spawn, it emulates adding tasks from different threads or/and modules
    for(int task=0; task <_NUM_TASKS; ++task) {
        cilk_spawn [&task_duration,task] {
            task_duration[task] = time([&] { // Start time measurment
                    compute_task();
                });
        } ();
    }
    cilk_sync;
#else
    std::cerr << "Executing TBB baseline" << std::endl;
    tbb::task_arena arena;
    tbb::task_group group;
    // though this is not optimal to use for-loop+task_group, it emulates adding tasks from different threads or/and modules
    for(int task=0; task <_NUM_TASKS; ++task) {
        arena.execute([&]{
            group.run([&,task]{
                task_duration[task] = time([&] { // Start time measurment
                    compute_task();
                });
            });
        });
    }

    arena.execute([&]{ group.wait(); });
#endif
}

void execute_flat_arena(duration_array_t& task_duration, const std::shared_ptr<arenas::base_arena>& arena, const arenas::topology& topology) {
    const int tokens = topology.num_nodes();
    arena->execute([tokens,&task_duration]{
        int curr_task = 0;
        // Create pipeline
	    tbb::parallel_pipeline(tokens,
            // Generate parallel _NUM_TASKS
            tbb::make_filter<void,int>( tbb::filter::serial, [&curr_task](tbb::flow_control& fc)-> int{
                if( curr_task<_NUM_TASKS ) {
                    return curr_task++;
                 } else {
                    fc.stop();
                    return 0;
             }} ) &
             tbb::make_filter<int,void>( tbb::filter::parallel, [&task_duration](int task) {
                 task_duration[task] = time([&] {
                    compute_task();
                 });
             })
        );
	});
}

void execute_hierarchical_arena(duration_array_t& task_duration, const std::shared_ptr<arenas::base_arena>& arena) {
    for(int task=0; task<_NUM_TASKS; ++task) {
        arena->enqueue([&,task]() {
            task_duration[task] = time([&] { // Start time measurment
                compute_task();
            }); // end of measure
        }); // end of enqueue
    }

	// Wait for task completion
    arena->wait();
}

int main(int argc, char* argv[], char* arge[])
{
#ifdef __TASK_BASELINE__
    data_vector_t src_vect(DATA_SIZE);
    data_vector_t dst_vect(DATA_SIZE);

    double* src = &src_vect[0];
    double* dst = &dst_vect[0];

    std::chrono::duration<double> base_task_time = time([&] {
    for(int i=0; i<DATA_SIZE; i+= BLOCK_SIZE)
    {
        int begin = i, end = std::min(i+BLOCK_SIZE, DATA_SIZE);
        #pragma noinline
        fill_input(src, begin, end);
    }
    for(int t=0; t<ITERATIONS; ++t)
    {
        for(int i=0; i<DATA_SIZE; i+=BLOCK_SIZE)
        {
            int begin = i, end = std::min(i+BLOCK_SIZE, DATA_SIZE);
            #pragma noinline
            compute_sqrt_block(src, dst, begin, end, DATA_SIZE);
        }
    }
    });

    std::cout << "base_task_time: " << base_task_time.count() << " [sec]";
    return 0;
#endif
    __itt_pause();

#ifdef __MIC__
    // Setup process affinity mask not to run on the system core
    cpu_set_t affinity;
    CPU_ZERO(&affinity);
    sched_getaffinity(0, sizeof(affinity), &affinity);
    // Check if MIC system core in use.
    // System core should be not used during benchmarking
    if ( CPU_ISSET(0, &affinity) ) {
        int cpu_count = CPU_COUNT(&affinity);
        std::cout << "setting application mask, cpu count = " << cpu_count << std::endl;
        CPU_ZERO(&affinity);
        for(int i=1;i<cpu_count-3;++i)
            //if(i&1)
                CPU_SET(i, &affinity);
        sched_setaffinity(0, sizeof(affinity), &affinity);
    }

    arenas::topology                    topology(MIC_TEAMS, 4);
#else
    arenas::topology                    topology;
#endif

    duration_array_t task_duration(_NUM_TASKS);

#if !(defined(__OMP_BASELINE__) || defined(__TBB_BASELINE__) || defined(__CILK_BASELINE__))
    std::shared_ptr<arenas::base_arena> arena;
    bool isFlat = (argc > 1) && strcmp(argv[1], "flat") == 0;
	if ( isFlat ) {
        std::cerr << "Executing with flat arena" << std::endl;
        arena.reset( new arenas::flat_arena(topology) );
    }
	else {
        std::cerr << "Executing with hierarchical arena" << std::endl;
        arena.reset( new arenas::hierarchical_arena(topology) );
	}
#else
    // Creat all worker threads
    initialize_threading();
#endif
    __itt_resume();

    std::chrono::duration<double> application_time = time([&]{
#if !(defined(__OMP_BASELINE__) || defined(__TBB_BASELINE__) || defined(__CILK_BASELINE__))
        if ( isFlat )
            execute_flat_arena(task_duration, arena, topology);
        else
            execute_hierarchical_arena(task_duration, arena);
#else
        execute_baseline(task_duration);
#endif
    }); // end of application measure

    std::sort(task_duration.begin(), task_duration.end());
    std::chrono::duration<double> total_time = std::accumulate(task_duration.begin(), task_duration.end(), std::chrono::duration<double>(0.));
    std::cout << "application_time: " << application_time.count() << " [sec] "
              << "total_time: " << total_time.count() << " [sec] "
              << "minimum_time: " << task_duration[0].count() << " [sec] "
              << "median_time: " << task_duration[_NUM_TASKS/2+1].count() << " [sec] "
              << "mean_time: " << total_time.count()/_NUM_TASKS << " [sec] "
              << "maximum_time: " << task_duration[_NUM_TASKS-1].count() << " [sec]" << std::endl;

    return 0;
}

