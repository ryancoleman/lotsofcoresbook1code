#!/bin/bash

runs=6
phi="7xxx"
mpss="3.4"
compiler="intel15.0.0"
mpi="intelMPI5.0.1.035"
postfix="xxx"

echo "# package size[byte], ranks, transfer time per package [us], transfer rate [Gbyte/s]" > datatransfer_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${mpi}_${postfix}.log

for i in 1 16
do
    len=1
    for l in `seq 1 1 25`
    do
	mpiicpc -openmp -DLENGTH=${len} -std=c++11 -DBATCH_OUTPUT -lrt -o dataTransfer.x dataTransfer.cpp
	mpirun -np ${i} ./dataTransfer.x > temp.log
	cat temp.log
	IFS=' ' read -a array <<< `cat temp.log`
	a=${array[0]}
	b=${array[1]}
	c=${array[2]}
	d=${array[3]}
	
	for j in `seq 1 1 $runs`
	do 
	    sleep 2
	    mpirun -np ${i} ./dataTransfer.x > temp.log
	    cat temp.log
    	    IFS=' ' read -a array <<< `cat temp.log`
	    if (( $(bc <<< "${array[2]} < $c") == 1 ))
	    then
		a=${array[0]}
		b=${array[1]}
		c=${array[2]}
		d=${array[3]}
	    fi
	done
	
	echo "$a  $b  $c  $d" >> datatransfer_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${mpi}_${postfix}.log
	echo "minimum: $a  $b  $c  $d"
	echo " "

	rm -f dataTransfer.x
	rm -f temp.log

	len=$(($len*2))
    done
    echo " " >> datatransfer_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${mpi}_${postfix}.log
    echo " " >> datatransfer_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${mpi}_${postfix}.log
done
