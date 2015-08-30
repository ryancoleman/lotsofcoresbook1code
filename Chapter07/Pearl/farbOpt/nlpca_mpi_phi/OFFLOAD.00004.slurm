#!/bin/bash
source ./common.sh
#!/bin/bash
#SBATCH -J offloadMPI        # Job Name
#SBATCH -o offloadMPI.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -N 4            # Total nodes
#SBATCH -n 4           # Total number of mpi tasks requested
#SBATCH -p normal-mic  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 00:10:00     # Run time (hh:mm:ss) - 1.5 hours

ibrun ./mpiTrain_$APP.off NameNotUsed.dat $APP.param mpiOffTiming
