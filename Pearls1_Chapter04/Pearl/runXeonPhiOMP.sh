#Script requires number of threads passed in as argument
echo "OMP_NUM_THREADS = $1"
mem=$(( (3500/$1)))
echo "mem == ${mem}m"
scp ./main.Intel.omp.exe mic0:~
scp ./inputs_SMC mic0:~
#may need to change compiler path depending on your setup
compilerPath="/opt/intel/composer_xe_2015"
ssh ${HOSTNAME}-mic0 "source $compilerPath/bin/compilervars.sh intel64; ulimit -s unlimited; export KMP_PLACE_THREADS=60c,4t; export KMP_STACKSIZE=${mem}m; export OMP_NUM_THREADS=$1;  export KMP_AFFINITY=balanced; export LD_LIBRARY_PATH=$compilerPath/lib/mic:$LD_LIBRARY_PATH; ~/main.Intel.omp.exe"

