#!/bin/bash

. /opt/intel/composerxe/bin/compilervars.sh intel64
export PATH=/opt/gcc-4.7.0/bin:${PATH}
export LD_LIBRARY_PATH=/opt/gcc-4.7.0/lib64:${LD_LIBRARY_PATH}
export SINK_LD_LIBRARY_PATH=$MIC_LD_LIBRARY_PATH
make clean

if [ "$1" == "knc" ]; then
	export cmd="micnativeloadex ./gems"
	export OUT_FILE=baseline_knc.csv
        export plat=knc

else
	export cmd="./gems"
	export OUT_FILE=baseline_xeon.csv
fi

touch build.cmd
make gems BASELINE=TASK PLAT=${plat} TASKS=1
task_base=`${cmd} | grep sec | cut -d ' ' -f 2`
echo "Task_Baseline,$task_base" > $OUT_FILE

echo "Baseline,Tasks,Application,Total,Minimum,Median,Average,Maximum" >> $OUT_FILE
for base in OPENMP TBB CILK
do
	for task in 1 2 3 4 5 10 15 20 30
	do
		touch build.cmd
		make gems BASELINE=${base} PLAT=${plat} TASKS=${task}
		res=`${cmd} | grep sec`
		echo $res
		apl=`echo $res | cut -d ' ' -f 2`
		tot=`echo $res | cut -d ' ' -f 5`
		min=`echo $res | cut -d ' ' -f 8`
		med=`echo $res | cut -d ' ' -f 11`
		avg=`echo $res | cut -d ' ' -f 14`
		max=`echo $res | cut -d ' ' -f 17`

		echo "$base,$task,$apl,$tot,$min,$med,$avg,$max" >> $OUT_FILE
	done
done

cat $OUT_FILE
