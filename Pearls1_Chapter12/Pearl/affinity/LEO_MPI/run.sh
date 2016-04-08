#!/bin/bash

groups=4
groupSize=4

if [[ "$#" == "2" ]]
then 
    groups=$1
    groupSize=$2
fi

make GROUP_SIZE=$groupSize

mpiString=""

for i in `seq 0 1 $(($groups-1))`
do
    if [[ "$i" == "0" ]]
    then
	mpiString="-np 1 -env MIC_KMP_PLACE_THREADS=$(($groupSize/4))c,${groupSize}t,$(($i*($groupSize/4)))O ./printThreadCoreAssignment.x"
    else
	mpiString="$mpiString : -np 1 -env MIC_KMP_PLACE_THREADS=$(($groupSize/4))c,${groupSize}t,$(($i*($groupSize/4)))O ./printThreadCoreAssignment.x"
    fi
done

# mic-env prefix
export MIC_ENV_PREFIX=MIC

mpirun $mpiString

make clean
