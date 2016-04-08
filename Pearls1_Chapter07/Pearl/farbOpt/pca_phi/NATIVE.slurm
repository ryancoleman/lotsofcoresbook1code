#!/bin/bash
#SBATCH -J nativePerf        # Job Name
#SBATCH -o nativePerf.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -N 1           # Total number of mpi tasks requested
#SBATCH -n 1           # Total number of mpi tasks requested
#SBATCH -p normal-mic  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 00:10:00     # Run time (hh:mm:ss) - 1.5 hours

sh RUN_NATIVE
