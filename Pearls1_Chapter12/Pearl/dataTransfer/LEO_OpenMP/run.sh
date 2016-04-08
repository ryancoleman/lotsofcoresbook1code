#!/bin/bash

runs=6
phi="7xxx"
mpss="3.4"
compiler="intel15.0.0"
postfix="fw"

echo "# package size[byte], ranks, transfer time per package [us], transfer rate [Gbyte/s]" > datatransfer_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${postfix}.log

for i in 1 16
do
    len=1
    for l in `seq 1 1 25`
    do
	icpc -openmp -DLENGTH=${len} -DNUM_THREADS=${i} -std=c++11 -DBATCH_OUTPUT -lrt -o dataTransfer.x dataTransfer.cpp
	./dataTransfer.x > temp.log
	cat temp.log
	IFS=' ' read -a array <<< `cat temp.log`
	a=${array[0]}
	b=${array[1]}
	c=${array[2]}
	d=${array[3]}
	
	for j in `seq 1 1 $runs`
	do 
	    sleep 2
	    ./dataTransfer.x > temp.log
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
	
	echo "$a  $b  $c  $d" >> datatransfer_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${postfix}.log
	echo "minimum: $a  $b  $c  $d"
	echo " "

	rm -f dataTransfer.x
	rm -f temp.log

	len=$(($len*2))
    done
    echo " " >> datatransfer_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${postfix}.log
    echo " " >> datatransfer_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${postfix}.log
done
