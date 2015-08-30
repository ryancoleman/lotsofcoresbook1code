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

#include <algorithm>
#include <random>
#include <functional>
#include <vector>
#include <assert.h>

#define TBB_PREVIEW_TASK_ARENA 1
#include <tbb/task.h>
#include <tbb/task_arena.h>
#include <tbb/atomic.h>
#include <tbb/partitioner.h>
#include <tbb/task_group.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/tbb_thread.h>
#undef min

#if TBB_INTERFACE_VERSION < 7005
#define current_thread_index current_slot
#endif

#define MIC_THREADS_PER_CORE 4  // Don't change, it's Xeon Phi(tm) hardware property

#include "observers.h"
#ifndef _WIN32
#ifndef __MIC__
#include <numa.h> // NUMA library is not supported on Intel(r) Xeon Phi(tm) platform
#endif
#else
#include <Systemtopologyapi.h>
#endif

namespace arenas {

class topology {

    typedef std::vector<std::vector<int>> machine_topology_t;

    machine_topology_t my_topology;
    std::vector<int>   my_logical_cores;

    void create_id_array(std::vector<int>& cpud_ids, int max_threads) {
#ifdef _WIN32
        DWORD_PTR process_mask = 0, system_mask;
        GetProcessAffinityMask(GetCurrentProcess(), &process_mask, &system_mask);
        for(int i=0; i<sizeof(KAFFINITY)*8 && cpud_ids.size()<max_threads; ++i) {
            if ( (((__int64)1<<i) & process_mask) != 0)
            {
                cpud_ids.push_back(i);
            }
        }
#else
        cpu_set_t process_mask;
        CPU_ZERO(&process_mask);
        sched_getaffinity(0, sizeof(process_mask), &process_mask);
        for(int i=0; i<CPU_SETSIZE && cpud_ids.size()<max_threads; ++i) {
            if ( CPU_ISSET(i, &process_mask) ) {
                cpud_ids.push_back(i);
            }
        }
#endif
   }

public:

    const machine_topology_t& operator()( void ) const { return my_topology;}
    const std::vector<int>& operator()(int i) const { return my_topology[i]; }
    const int& operator[](int i) const { return my_logical_cores[i]; }

    int num_cpus() const  { return (int)my_logical_cores.size(); }
    int num_nodes() const { return (int)my_topology.size(); }

    topology() {
#ifdef _WIN32
        ULONG num_nodes = 0;
        BOOL res = GetNumaHighestNodeNumber(&num_nodes);
        my_topology.resize(++num_nodes);
        for(USHORT n = 0; (n < num_nodes) && res; n++) {
            GROUP_AFFINITY node_affinity;
            res = GetNumaNodeProcessorMaskEx(n, &node_affinity);
            if ( !res ) break;
            for(int i=0;i<sizeof(KAFFINITY)*8;++i) {
                if ( (((__int64)1<<i) & node_affinity.Mask) != 0)
                {
                    my_topology[n].push_back(i);
                    my_logical_cores.push_back(i);
                }
            }
        }

#elif !defined(__MIC__)
        const int max_cpus = numa_num_possible_cpus();
        struct bitmask *mask = numa_allocate_cpumask();
        const int num_nodes = numa_num_task_nodes();

        // Check against process affinity as well
        cpu_set_t process_mask;
        CPU_ZERO(&process_mask);
        sched_getaffinity(0, sizeof(process_mask), &process_mask);
     
        std::cout << "CPU count: " << numa_num_task_cpus() << "  NUMA nodes: " << num_nodes << std::endl;

        my_topology.resize(num_nodes);
        for(int n = 0; n < num_nodes; n++) {
            int err = numa_node_to_cpus(n, mask);
            if( err ) { perror("numa_node_to_cpus"); exit( EXIT_FAILURE ); }
            std::cout << "\tNode #" << n << " CPU ids:";
            for(int i=0; i<max_cpus; i++)
                if ( numa_bitmask_isbitset(mask, i) && CPU_ISSET(i, &process_mask) ) {
                    my_topology[n].push_back(i);
                    my_logical_cores.push_back(i);
                    std::cout << ' ' << i;
                }
            std::cout << std::endl;
        }
        numa_free_cpumask(mask);
#endif

        // If NUMA API fails, request node information from
        // the process affinity mask
        if ( my_topology.empty() || my_topology[0].empty() ) {
            create_id_array(my_logical_cores, INT_MAX);
            my_topology.resize(1);
            my_topology[0] = my_logical_cores;
        }
    }

    topology(int max_threads) {
        create_id_array(my_logical_cores, max_threads);
        my_topology.push_back(my_logical_cores);
    }

    // Parameter: number of threads per core to be used by a team
    topology(int num_teams, int thread_core_team) {
        create_id_array(my_logical_cores, INT_MAX);
        my_topology.resize(num_teams);
        const int cores_per_team = (my_logical_cores.size()/MIC_THREADS_PER_CORE)/num_teams;
        const int thread_per_team = thread_core_team * cores_per_team;

        int num_threads = num_teams*thread_per_team;

        for(int team=0;team<num_teams;++team)
        {
            int core_offset = team*cores_per_team;
            for(int thread=0;thread<thread_per_team;++thread)
            {
                int core_id = core_offset + thread / thread_core_team;
                int thread_id = thread % thread_core_team;

                int thread_num = core_id*MIC_THREADS_PER_CORE+thread_id; // thread number in over all threads
                my_topology[team].push_back( my_logical_cores[thread_num]);
            }
        }
    }

    ~topology() {
    }
};

/*! /brief Submission arena base interface
 * Arena represents a place where in a work/tasks are introduced.
 * Worker threads join an arena in order to execute introduced work.
 * Each arena has its own concurrency level which might be different
 * from the default value.
 * Current implementation supports:
 *    1. single level or flat arena - intended for usage in small scale systems or
 *             when NUMA related effects are not important
 *    2. dual level or hierarchical arena - intended for usage in large scale systems
 *             such as Intel(r) Xeon Phi(tm), or when NUMA related effects are important                         
*/
class base_arena {
public:
    virtual ~base_arena() {
    }

    //! Wait for all work to be completed in an arena
    void wait() {
        my_arena.execute([&]{ my_group.wait(); } );
    }

    //! Submit a work represented by a lambda function in to an arena.
    //! A representative task is enqueued for later execution
    /*!
      \param functor A lambda function to be submitted.
      \param teamId  A team number to which the work should be enqueued. Relevant only for hierarchical arena
    */
    template<typename F>
    void enqueue(const F& _f, int teamId = -1) {
        enqueue_int(_f, teamId);
    }

    //! Execute a work represented by a lambda function in an arena.
    //! A representative task is executed immediately, calling thread joins execution
    /*!
      \param functor A lambda function to be executed.
      \param teamId  A team number on which the work should be executed. Relevant only for hierarchical arena
    */
    template<typename F>
    void execute(const F& _f, int teamId = -1) {
        execute_int(_f, teamId);
    }

protected:
#ifdef __ENABLE_TRAPPING__
    // DISCLAIMER: this technique is not recommended by TBB and should be avoided in production
    class arena_trapper
    {
    protected:
        struct arena_traper_task : public tbb::task {
            arena_trapper& owner;
            int numThreads;
            arena_traper_task( arena_trapper& _owner, int n ) : owner(_owner), numThreads(n) {}

            tbb::task* execute () {
                owner.barrier++;
                while ( owner.barrier < numThreads )
                    tbb::this_tbb_thread::yield(); // Wait until all workers are ready

                // Bundle the worker thread to arena, wait for incoming work
                // Master notifies when workers should leave
                owner.root->wait_for_all();

                // Wait until all workers 
                owner.barrier--;
                while ( owner.barrier > 0 )
                    tbb::this_tbb_thread::yield();
                return NULL;
            }
        };

    public:
        arena_trapper()
        : wait_context(tbb::task_group_context::isolated, tbb::task_group_context::default_traits
                       | tbb::task_group_context::concurrent_wait)
        , root(NULL)
        {}

        //! Execute a work represented by a lambda function in an arena.
        //! A representative task is executed immediately, calling thread joins execution
        /*!
            \param _cpu_set A vector of logical core ids where to locate worker threads
         */
        void init(int numThreads) {
            root = new ( tbb::task::allocate_root(wait_context) ) tbb::empty_task;
            root->set_ref_count(2);
            barrier = 0;
            for ( int i = 1; i < numThreads; ++i )
                tbb::task::spawn( *new(tbb::task::allocate_root(wait_context))
                                      arena_traper_task(*this, numThreads) );
            barrier++;
            while ( barrier < numThreads )
                tbb::this_tbb_thread::yield();
        }

        void release() {
            if ( barrier == 0 ) return;
            root->decrement_ref_count();
            barrier--;
            while ( barrier > 0 )
                    tbb::this_tbb_thread::yield();
            tbb::task::destroy(*root);
        }

    protected:
        tbb::task_group_context wait_context;
        tbb::task*              root;
        tbb::atomic<long>       barrier;
    };

    arena_trapper  trapper;
#endif

    typedef std::function<void(void)> arena_functor;

    base_arena() : my_arena(), my_group() {}

    virtual void enqueue_int(const arena_functor& func, int teamId) = 0;
    virtual void execute_int(const arena_functor& func, int teamId) = 0;

    tbb::task_arena                 my_arena;
    tbb::task_group                 my_group;
};

/*! /brief A class implements a flat, one level arena
 */
class flat_arena : public base_arena {
    observers::pinning_observer     my_observer;

public:
    flat_arena(const topology& _t)
    : my_observer(my_arena, &_t[0], true)
    {
        int num_threads = _t.num_cpus();

        // Initialize the arena class to desired concurrency,
        // leave one arenas slot for the master thread
        my_arena.initialize(num_threads, 1);

#ifdef __ENABLE_AFFINITY__
        my_observer.observe(true);
#endif

#ifdef __ENABLE_TRAPPING__
        // Initialized threads in top level arena
        my_arena.execute([&]{
            trapper.init(num_threads);
        });
#endif

#if 0  // !!!! Remove this code
        my_arena.execute([&]{
            // Worker initialization, warm-up threads and set affinity of workers
            tbb::atomic<int> worker_count;
            worker_count = num_threads;
            const int* cpu_id = &observers::topology::get().cpu_ids[0];
            tbb::parallel_for(0, num_threads, 1, [&](int i) {
                cpu_set_t affinity;
                CPU_ZERO(&affinity);
                CPU_SET(cpu_id[i], &affinity);
                std::cout << "set affinity to:" << cpu_id[i] << std::endl;

                sched_setaffinity(0, sizeof(affinity), &affinity);
                // DISCLAIMER: this barrier is against TBB guarantees and should not be used in production
                worker_count--;
                while(worker_count > 0)
                    tbb::this_tbb_thread::yield();
            }, tbb::simple_partitioner() );
        });
#endif
}

   ~flat_arena() {
#ifdef __ENABLE_TRAPPING__
            trapper.release();
#endif
   }

protected:
    // Inherited methods
    void enqueue_int(const arena_functor& func, int teamId)
    {
        my_arena.execute([&] {
            my_group.run(func);
        });
    }

    void execute_int(const arena_functor& func, int teamId)
    {
        my_arena.execute(func);
    }

};

/*! /brief A class implements a dual level hierarchical arena
 */
class hierarchical_arena : public base_arena {
public:
    hierarchical_arena(const topology& _t)
    :   arena_array(NULL)
    ,   master_mask(_t.num_nodes())
    ,   my_observer(my_arena, &master_mask[0], true)
    {
        arena_count = master_mask.size();
        my_arena.initialize(arena_count, 1);  // Allocate 1 slot for each "execution" arena, enable master thread to join

        arena_array = new nested_arena[arena_count];
        for(int i=0; i<arena_count; ++i) {
            master_mask[i] = _t(i)[0];
        }

#ifdef __ENABLE_AFFINITY__
        my_observer.observe(true);
#endif

#ifdef __ENABLE_TRAPPING__
        // Initialized threads in top level arena
        my_arena.execute([&]{
            trapper.init(arena_count);
        });
#endif

#if 0 // It looks like we should not use barrier at all
        // Initialize second level arenas
        tbb::atomic<long> barrier;
        barrier = arena_count;
#endif
        my_arena.execute([&]{
            tbb::parallel_for(0, arena_count, [&](int i)
            {
                int slot = i;
                arena_array[slot].init(_t(i));
#if 0
                barrier--;
                // DISCLAIMER: this barrier is against TBB guarantees and should not be used in production
                while(barrier > 0)
                    tbb::this_tbb_thread::yield();
#endif
            }, tbb::simple_partitioner());
        });        
    }

    ~hierarchical_arena()
    {
        delete[] arena_array;
    }

protected:

    class nested_arena : public tbb::task_arena {
    public:
        nested_arena() : my_observer(*this, 0, false) {}
        ~nested_arena() { release(); }

        void init(const std::vector<int>& _cpu_set)
        {
            size_t set_size = _cpu_set.size();
            // Initialize arena with desired number of threads, leave one slot for master
            tbb::task_arena::initialize((int)set_size, 1);
#ifdef __ENABLE_AFFINITY__
            my_observer.set_map(&_cpu_set[0]);
            my_observer.observe(true);
#endif

#ifdef __ENABLE_TRAPPING__
            tbb::task_arena::execute([&]{
                trapper.init(set_size);
            });
#endif
        }

        void release() {
#ifdef __ENABLE_TRAPPING__
            trapper.release();
#endif
        }

    protected:
#ifdef __ENABLE_TRAPPING__
        arena_trapper trapper;
#endif
        observers::pinning_observer my_observer;
    };

    void enqueue_int(const arena_functor& _func, int teamId) {
        // Need to allocate slot wherein to put a new task
        my_arena.execute([&] {
            my_group.run( [&,_func, teamId] {
                this->execute_int(_func, teamId);
            });
        });
    }

    // inherited functions implementation
    void execute_int(const arena_functor& _func, int teamId)
    {
        if ( -1 == teamId ) {
            // Discover specific arena to execute on
            int arena_id = tbb::task_arena::current_thread_index();
            arena_array[arena_id].execute(_func);
        }
        else if ( -2 == teamId ) {
            my_arena.execute(_func);
        }
        else {
            assert( teamId < arena_count && "" );
            arena_array[teamId].execute(_func);
        }
    }

protected:
    nested_arena*  arena_array;
    int            arena_count;
    std::vector<int> master_mask;
    observers::pinning_observer my_observer;
};
} // namespace arenas
