========================================================================
    Gems_Code - build and execution instructions
========================================================================

Pre-requisites
--------------
1. Intel(R) Composer XE 2013 SP1 update 3 or later
2. GCC 4.7 (later GCC versions require Intel(R) Compser XE 2015)
3. NUMA development package (libnuma-dev)
4. Intel(R) Many Integrated Core (MIC) Platform Software Stack (MPSS)
   version 3.2.3 or later
5. Bash shell

Build and Makefile parameters
-----------------------------
Make utility file provided within a project is used to build the benchmark
for Intel(r) Xeon(r) and Intel(r) Xeon Phi(tm) hardware. In addition,
default invocation (w/o make parameters) invokes executables.

Makefile parameters:
- BASELINE= specifies the threading library for baseline benchmarking,
            supported values: OPENMP, CILK, TBB. If value is not specified,
            arena's based code is executed.

- PLAT=     select platform to build for. Default value is 'avx',
            build for Intel(r) Xeon(r). 'knc' - build for Intel(r)
            Xeon Phi(tm)
- TASKS=   specifies the number of tasks to execute. For AVX platform
           default value is 60, for KNC platform the default is 20.
- TEAMS=   number of live tokens / hierarchical arena groups to be
           created when running with arenas. Default value 20.

Benchmark execution
-------------------
Benchmark scripts:
- run_mic.sh, this script executes performance evaluation of arena's
              related benchmarks on Intel(r) Xeon Phi(tm). It invokes
              the build process for 1,2,3,4,5,10,15,20,30 and 60 team,
              while TASKS=120.
              When 'flat' is specified as an argument, the flat arena
              approach is executed; otherwise, hierarchical arena is used.
              The result are store as coma separated values (csv) in
              flat_run.csv and hierarchical_run.csv respectively.
- run_baseline.sh, is used for execution of the baseline measurements
               for the threading libraries. When 'knc' is specified as
               an argument the benchmarks is invoked on Intel(r) Xeon
               Phi(tm). The results of Intel(r) Xeon Phi(tm) are stored
               in run_baseline_knc.csv, and for Intel(r) Xeon(r)
               invocation results are stored in run_baseline.csv

NUMA evaluation:
For NUMA evaluation of hierarchical arenas the following commands
shall be executed:

make clean
make gems TEAMS=2
./gems

/////////////////////////////////////////////////////////////////////////////
