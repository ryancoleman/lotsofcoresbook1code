#!/bin/bash
source ./common.sh
#!/bin/bash
#SBATCH -J hostMPI        # Job Name
#SBATCH -o hostMPI.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -N 2            # Total nodes
#SBATCH -n 2           # Total number of mpi tasks requested
#SBATCH -p normal  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 00:10:00     # Run time (hh:mm:ss) - 1.5 hours

ibrun ./mpiTrain_$APP.omp NameNotUsed.dat $APP.param mpiHostTiming
