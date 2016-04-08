#Script requires number of threads passed in as argument
ulimit -s unlimited
echo "OMP_NUM_THREADS = $1"
mem=$(( (3500/$1)))
echo "mem == ${mem}m"
#may need to change compiler path depending on your setup
compilerPath="/opt/intel/composer_xe_2015"
mpiexec.hydra -np 1 -env KMP_STACKSIZE=${mem}m -env OMP_NUM_THREADS=$1 -env KMP_AFFINITY=scatter -env LD_LIBRARY_PATH=$compilerPath/lib:$LD_LIBRARY_PATH ./main.Intel.mpi.omp.exe 

