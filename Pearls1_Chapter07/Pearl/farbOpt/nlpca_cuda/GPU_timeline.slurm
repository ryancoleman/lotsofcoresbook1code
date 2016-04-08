#!/bin/bash
#SBATCH -J gpuTimeline        # Job Name
#SBATCH -o gpuTimeline.o%j    # Output and error file name (%j expands to jobID)
#SBATCH -N 1           # Total number of mpi tasks requested
#SBATCH -n 1           # Total number of mpi tasks requested
#SBATCH -p gpu  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 00:10:00     # Run time (hh:mm:ss) - 1.5 hours

. ./common.sh

./gen_$APP $APP.train.dat $EXAMPLES $VAR 1234

nvprof --output-profile gpu.timeline.nvprof ./train_$APP.cuda $APP.train.dat $APP.param

./gen_$APP $APP.pred.dat 1000 0 1
./pred_$APP $APP.param $APP.pred.dat > output.txt

# create file for gnuplot
tail -n +3 output.txt > plot.txt

#rm *.dat
