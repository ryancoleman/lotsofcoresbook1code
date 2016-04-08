#!/bin/bash

. /opt/intel/composerxe/bin/compilervars.sh intel64
export PATH=/opt/gcc-4.7.0/bin:${PATH}
export LD_LIBRARY_PATH=/opt/gcc-4.7.0/lib64:${LD_LIBRARY_PATH}
export SINK_LD_LIBRARY_PATH=$MIC_LD_LIBRARY_PATH

make clean

if [ "$1" == "flat" ]; then
	export FLAT_PARAM="-a flat"
	export OUT_FILE=flat_run.csv
else
	export OUT_FILE=hierarchical_run.csv
fi

echo "Teams,Application,Total,Minimum,Median,Average,Maximum" > $OUT_FILE
for teams in 1 2 3 4 5 10 15 20 25 30 
do
	touch build.cmd
	make gems PLAT=knc TEAMS=$teams TASKS=120
	res=`micnativeloadex gems ${FLAT_PARAM} | grep sec`
	echo $res
	apl=`echo $res | cut -d ' ' -f 2`
	tot=`echo $res | cut -d ' ' -f 5`
	min=`echo $res | cut -d ' ' -f 8`
	med=`echo $res | cut -d ' ' -f 11`
	avg=`echo $res | cut -d ' ' -f 14`
	max=`echo $res | cut -d ' ' -f 17`

	echo "$teams,$apl,$tot,$min,$med,$avg,$max" >> $OUT_FILE
done

echo "Tasks,Application,Total,Minimum,Median,Average,Maximum" >> $OUT_FILE
for tasks in 1 5 10 15 20 30 40 50 60 70 80 90 100 110 120
do
        touch build.cmd
        make gems PLAT=knc TEAMS=20 TASKS=$tasks
        res=`micnativeloadex gems ${FLAT_PARAM} | grep sec`
        echo $res
        apl=`echo $res | cut -d ' ' -f 2`
        tot=`echo $res | cut -d ' ' -f 5`
        min=`echo $res | cut -d ' ' -f 8`
        med=`echo $res | cut -d ' ' -f 11`
        avg=`echo $res | cut -d ' ' -f 14`
        max=`echo $res | cut -d ' ' -f 17`

        echo "$tasks,$apl,$tot,$min,$med,$avg,$max" >> $OUT_FILE
done

cat $OUT_FILE
