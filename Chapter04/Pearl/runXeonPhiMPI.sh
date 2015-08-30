#Script requires number of threads passed in as argument
echo "OMP_NUM_THREADS = $1"
mem=$(( (3500/$1)))
echo "mem == ${mem}m"
scp ./main.Intel.mpi.omp.exe mic0:~
scp ./inputs_SMC mic0:~
#may need to change compiler and MPI path depending on your setup
compilerPath="/opt/intel/composer_xe_2015"
MPIPath="/opt/intel/impi/5.0.1.035/mic"
ssh ${HOSTNAME}-mic0 "source $compilerPath/bin/compilervars.sh intel64; source $MPIPath/bin/mpivars.sh; ulimit -s unlimited; mpiexec.hydra -np 1 -env KMP_PLACE_THREADS=60c,4t -env KMP_STACKSIZE=${mem}m -env OMP_NUM_THREADS=$1 -env KMP_AFFINITY=balanced -env LD_LIBRARY_PATH=$compilerPath/lib/mic:$LD_LIBRARY_PATH ~/main.Intel.mpi.omp.exe"

