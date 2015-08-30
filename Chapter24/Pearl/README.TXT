The files in these directory supplement the article 
"Profiling-Guided Optimization" (2014, Andrey Vladimirov) 
published as a chapter in the book
"High Performance Parallelism Pearls – Successful Approaches 
for Multicore and Many-core Programming", 2014.

QUICK START:
                 make
                 make run-cpu
                 make run-mic

IMPORTANT: prior to running the MIC executable, the OpenMP library
           for the MIC architecture must be copied to, or shared with 
           the coprocessor. One of the ways to do it is by running:
           sudo scp /opt/intel/composerxe/lib/mic/libiomp5.so mic0:/lib64/

DETAILS:

Directories 1/, 2/, 3/, 4/ and 5/ contain the source code and make files
of the code at various levels of optimization, as discussed in the article.

1/ contains serial transposition code
2/ contains transposition code parallelized using OpenMP
3/ contains code optimized with loop tiling
4/ contains further optimized code in which the inner loop is regularized by partitioning 
   the matrix into regions
5/ contains the final result, with additional parallelism exposed by means of planning

Each of these directories contains four files:
  Transpose.cc - the code of the transposition function
  Transpose.h  - the header file for the corresponding .cc file
  Main.cc      - the driver that performs matrix initialization and benchmarks the code
  Makefile     - the makefile supporting target "all"

Running "make" in each directory produces two executables:
  runme-CPU    - executable for Intel Xeon architecture
  runme-MIC    - native executable for Intel Xeon Phi (MIC) architecture

In order to run the CPU executable, pass two command line arguments to it:
matrix size and number of trials (at least 10). For example:
               ./runme-CPU 4000 32

In order to run the MIC executable, copy it to the coprocessor using scp
and then run in a similar fashion to the CPU executable:
               scp runme-CPU mic0:~/
               ssh mic0 ./runme-CPU 4000 32

See Makefile in the current directory for optimized environment variable
settings (number of threads and thread affinity controls).
