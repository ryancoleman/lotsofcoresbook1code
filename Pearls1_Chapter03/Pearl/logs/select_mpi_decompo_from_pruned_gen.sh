#!/bin/bash 

# Note that we  might get several choices for a decomposition using $N task and
# the idea here is to take the best. Note that we might not get a complete set.
grep nprocj options.nml || exit
ls decompo_stats_*.txt || exit
N=`grep nprocj options.nml |awk '{print $NF}'`
FILE=`ls decompo_stats_*.txt`
echo $FILE 
echo $N
test -d NEW && ls -l NEW/mpi*|wc -l
test -d NEW && rm -fr NEW
mkdir -p NEW
while [ $N -ge 1 ]
do
  echo "Finding the best decomposition for $N tasks"
  echo "Debug print showing candidates"
  grep -C2 "non-empty tasks: *$N\$"  $FILE | tr  "\n\r" " "|sed 's;=====*;\n\r";g'|sed 's;" Auto-decompo: *;;g'|sed 's;1 non-empty tasks: ;;g'|sed 's;area.*;;g'|sed 's; global Cmax/Cmin:   ;;g' |grep '[0-9]'|awk '{ print $0| "sort -k 2" }'
  txt=`grep -C2 "non-empty tasks: *$N\$"  $FILE | tr  "\n\r" " "|sed 's;=====*;\n\r";g'|sed 's;" Auto-decompo: *;;g'|sed 's;1 non-empty tasks: ;;g'|sed 's;area.*;;g'|sed 's; global Cmax/Cmin:   ;;g' | tr "\r" " " |grep '[0-9]'|awk '{ print $0| "sort -k 2" }'|head -1`
#echo "txt is $txt"
  XXX=`echo $txt|awk '{print $2 " " $4}'`
  echo "WINNER $XXX"
  XX=`echo $txt|awk '{print $2 " " $5}'`
  echo "ratio $XX"
  ratio=`echo $XX|awk '{print $NF}'`
  from=`echo $txt|awk '{print $1}'`
  to=`echo $txt|awk '{print $2}'`
  echo "JWP: $ratio"
# FIXME hardcoding
  if (( $(bc <<< "$ratio <= 1.4") == 1 )); then 
    echo "YES, jwp - ratio <= 1.4: $ratio"
    echo "cp mpi_decompo_${from}x1.txt NEW/mpi_decompo_${to}x1.txt"
    cp mpi_decompo_${from}x1.txt NEW/mpi_decompo_${to}x1.txt
  fi
  N=$(( $N - 1 ))
done
test -d NEW && ls -l NEW/mpi*|wc -l
exit
