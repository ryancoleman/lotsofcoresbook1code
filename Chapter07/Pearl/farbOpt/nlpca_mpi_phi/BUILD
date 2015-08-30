. ./common.sh
NO_IO_TEST=" -DNO_IO_EXAMPLES=$EXAMPLES "

FLAGS="-D$ACTIVATION -std=c99 -O3 -openmp -DMPI_NUM_COPROC_PER_NODE=1 "
FLAGS+=$MAX_RUN
FLAGS+=$NO_IO_TEST

INC=" -I $NLOPT_DIR/include  -I $MPICH_HOME/include "
LIB=" -L $NLOPT_DIR/lib -L $MPICH_HOME/lib "
#INC=" -I $NLOPT_DIR/include  -I /usr/include/mpich"
#LIB=" -L $NLOPT_DIR/lib -L /usr/lib/mpich "

icc $FLAGS genData.c -o gen_$APP

icc $FLAGS -Wno-unknown-pragmas -xhost ../mpiTrain.c -I . $INC $LIB -lmpich -lnlopt -lm -o mpiTrain_$APP.off

icc $FLAGS -Wno-unknown-pragmas -no-offload ../mpiTrain.c -I . $INC \
				$LIB -lnlopt -lmpich -lm -o mpiTrain_$APP.omp

icc $FLAGS -Wno-unknown-pragmas $NO_IO_TEST -no-offload ../pred.c -I . -lm -o pred_$APP

#FLAGS+=" -mmic -Wno-unknown-pragmas"
#NLOPT_DIR+="_native"

#icc $FLAGS ../mpiTrain.c -I . $INC $LIB -lnlopt -lmpich -lm -o mpiTrain_$APP.mic

