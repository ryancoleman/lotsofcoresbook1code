#!/bin/bash
#SBATCH -J hostPerf        # Job Name
#SBATCH -o hostPerf.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -N 1           # Total number of mpi tasks requested
#SBATCH -n 1           # Total number of mpi tasks requested
#SBATCH -p normal-mic  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 00:10:00     # Run time (hh:mm:ss) - 1.5 hours

sh RUN_HOST
