#!/bin/bash
#PBS -N LLT
#PBS -A your_account 
#PBS -j oe
#PBS -l nodes=4:ppn=2,walltime=02:00:00

#cd $PBS_O_WORKDIR

#performance tuning on MIC
export MIC_ENV_PREFIX=MIC
export MIC_USE_2MB_BUFFERS=2M
export MIC_KMP_AFFINITY=explicit,granularity=fine,proclist=[1-236:1] 
export KMP_AFFINITY=granularity=fine,compact,1,0

#initialize only the MIC to be used
export OFFLOAD_INIT=on_offload
#MIC list
export OFFLOAD_DEVICES=0,1

export PDLLT_N=5000
# MEMSIZE in unit of K wordsize
export PDLLT_MEMSIZE=2000
export PDLLT_NB=512
export PDLLT_P=2
export PDLLT_Q=1


#record these setting to output
echo 
#echo PBS_NUM_NODES=$PBS_NUM_NODES
#echo PBS_NUM_PPN=$PBS_NUM_PPN

echo MIC_ENV_PREFIX=$MIC_ENV_PREFIX
echo MIC_USE_2MB_BUFFERS=$MIC_USE_2MB_BUFFERS
echo MIC_KMP_AFFINITY=$MIC_KMP_AFFINITY
echo KMP_AFFINITY=$KMP_AFFINITY

echo OFFLOAD_INIT=$OFFLOAD_INIT
echo OFFLOAD_DEVICES=$OFFLOAD_DEVICES

echo PDLLT_N=$PDLLT_N
echo PDLLT_MEMSIZE=$PDLLT_MEMSIZE
echo PDLLT_NB=$PDLLT_NB
echo PDLLT_P=$PDLLT_P
echo PDLLT_Q=$PDLLT_Q

#check status of all mic (although we may not use all of them)
if [ $(mpirun -n $((PBS_NUM_NODES*PBS_NUM_PPN)) micctrl -w | grep online | wc -l) == $((PBS_NUM_NODES*PBS_NUM_PPN*2)) ]
then
echo mic check ok
mpirun -np $((PDLLT_P*PDLLT_Q)) ./pdlltdriver3.exe | grep "test result"  
mpirun -np $((PDLLT_P*PDLLT_Q)) ./pdlltdriver4.exe | grep "test result"
echo job done
else 
echo mic check failed
#exit -1
fi
date
#exit 0
