#Script requires number of threads passed in as argument
ulimit -s unlimited
echo "OMP_NUM_THREADS = $1" 
export OMP_NUM_THREADS=$1
mem=$(( (3500/$1)))
echo "mem == ${mem}m"
#may have to change LD_LIBRARY_PATH depending on your environment
export LD_LIBRARY_PATH=/opt/intel/composer_xe_2015/lib:$LD_LIBRARY_PATH
export KMP_STACKSIZE=${mem}m
export KMP_AFFINITY=scatter
./main.Intel.omp.exe
