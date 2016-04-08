#!/bin/bash
source ./common.sh
#!/bin/bash
#SBATCH -J gpuMPI        # Job Name
#SBATCH -o gpuMPI.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -N 4           # Total number of mpi tasks requested
#SBATCH -n 4           # Total number of mpi tasks requested
#SBATCH -p gpu  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 00:10:00     # Run time (hh:mm:ss) - 1.5 hours

ibrun ./mpiTrain_$APP.cuda NameDoesNotMatter.dat $APP.param gpuTiming
