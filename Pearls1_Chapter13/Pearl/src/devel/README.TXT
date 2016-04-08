=== Description ===
This directory contains the source code of step-by-step optimization of a basic
Zheng-Shu-Chew 3D LBM from a naive sequential implementation to an optimized 
OMP parallel implementation, and includes an offload version for execution in 
the Intel Xeon Phi.

In order of complexity the codes inside this directory are:

1. sequential  : naive sequential implementation
2. omp_v1      : added OMP to collision.F90
3. omp_v2      : added OMP to stream.F90, postcollision.F90, and poststream.F90
4. soa         : changed from array of sctructures to structure of arrays
5. soa_aligned : as soa version, but ensuring alignment to 64 byte boundary
6. offload     : offload version using both CPU and MIC coprocessor

These are development codes that eventually led to LBS3D and LBS3D-MPI. They 
are made available for educational purposes, and are not intended as production 
software. For more details see:

C. Rosales, Porting to the Intel Xeon Phi: Opportunities and Challenges. 
Extreme Scaling Workshop 2013 (XSCALE13), Boulder, CO (2013). 

2013-11-12 Carlos Rosales Fernandez
