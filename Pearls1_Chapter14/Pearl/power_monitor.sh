#!/bin/bash

while : ; do

    # Build the Timestamp string
    STR=`echo $(date +%N) | sed 's/^0*//'`
    STR=`printf "%03d" $(($STR/1000000))`
    STR=`echo $(date +%H:%M:%S)"."$STR`
    TIME="Timestamp: $STR"
    
    # Format the power and thermal data string as desired
    DATA="$({ micsmc -f; micsmc -t; sleep .1; } | grep -e Cpu -e temp -e Total -e freq | paste -d ' ' - -  | awk '{print $1,$4,"\t",$6,$7 }')"

   # Display our power and thermal data strings
    clear
    echo "$TIME"
    echo "$DATA"

done
