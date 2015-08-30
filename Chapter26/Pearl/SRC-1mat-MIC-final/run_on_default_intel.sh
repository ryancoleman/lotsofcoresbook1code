#!/bin/bash
source /opt/intel/composer_xe_2015/bin/compilervars.sh intel64
source /opt/intel/impi/5.0.1.035/bin64/mpivars.sh
make clean ; make pdlltdriver3.exe ; make pdlltdriver4.exe
. ./run_llt_mic.sh 
