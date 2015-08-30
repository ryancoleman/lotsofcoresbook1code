#!/bin/bash

nmics=0
host=$(hostname)

export I_MPI_MIC=1
export I_MPI_MIC_POSTFIX=.MIC
export I_MPI_FABRICS=tcp
source /opt/intel/impi_latest/intel64/bin/mpivars.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MIC_LD_LIBRARY_PATH
export KMP_AFFINITY=compact,granularity=fine
source /opt/intel/itac_latest/bin/itacvars.sh
export VT_LOGFILE_FORMAT=stfsingle

if [ "$1" == "cpu" ]; then

echo "$host:2" > machines
cp options ~/options

I_MPI_FABRICS=tcp \
VT_LOGFILE_FORMAT=stfsingle \
I_MPI_PIN_DOMAIN=omp \
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MIC_LD_LIBRARY_PATH} \
mpirun -trace \
    -machinefile machines -np 2 -env I_MPI_PIN=0 ~/options

cp options.single.stf options-cpu.single.stf

elif [ "$1" == "mic" ]; then

cp options ~/options
for imic in $(seq 0 $nmics)
do 
scp options.MIC mic${imic}:~/options.MIC
done

rm -rf machines
echo "$host" > machines # Boss rank will run on the host system
for imic in $(seq 0 $nmics)
do 
echo mic${imic} >> machines
done


I_MPI_FABRICS=tcp \
I_MPI_MIC=1 \
VT_LOGFILE_FORMAT=stfsingle \
I_MPI_PIN_DOMAIN=omp \
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MIC_LD_LIBRARY_PATH} \
mpirun -trace \
    -machinefile machines -env I_MPI_PIN=0 ~/options

cp options.single.stf options-mic.single.stf

elif [ "$1" == "both" ]; then

cp options ~/options
for imic in $(seq 0 $nmics)
do 
scp options.MIC mic${imic}:~/options.MIC
done

echo "$host:2" > machines
for imic in $(seq 0 $nmics)
do 
echo mic${imic} >> machines
done

I_MPI_FABRICS=tcp \
I_MPI_MIC=1 \
VT_LOGFILE_FORMAT=stfsingle \
I_MPI_PIN_DOMAIN=omp \
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MIC_LD_LIBRARY_PATH} \
mpirun -trace \
    -machinefile machines -env I_MPI_PIN=0 ~/options

cp options.single.stf options-both.single.stf

fi

echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~DONE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
