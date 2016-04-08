This version of miniMD was written as part of a chapter for 
"High Performance Parallelism Pearls".

It introduces a number of changes, compared to the original miniMD:
- Separate "local" and "remote" neighbor lists, allowing for the use of Newton's third law without atomics.
- The (rsq < cutforcesq) branch is replaced by a simple masked assignment to the force variable.
- A random "sort", to demonstrate the runtime improvement that results from sorting.
- Optimized versions of the force compute and neighbour list build steps, using intrinsics.
- A faster neighbor list build, that re-uses a transposed (AoS->SoA) version of a bin's stencil for all atoms in the bin.

To build the original code:
cp ORIGINAL/* .
make clean intel KNC=no ANSI_ALIAS=yes SP=yes
make clean intel KNC=yes ANSI_ALIAS=yes SP=yes

To build the optimized code with intrinsics (for AVX):
cp AVX/* .
make clean intel KNC=no ANSI_ALIAS=yes SP=yes

To build the optimized code with intrinsics (for KNC):
cp IMCI/* .
make clean intel KNC=yes ANSI_ALIAS=yes SP=yes

To run the code:
KMP_AFFINITY=compact,granularity=fine ./miniMD_intel -t [# threads] -gn 0 --half_neigh 1

Other combinations of precision, ghost newton and half neigh options are not supported by this version of miniMD.
Combinations of source files from the ORIGINAL/AVX/IMCI directories has not been tested and is not expected to work.
