APP=pca
FLAGS="-DUSE_LINEAR -std=c99 -O3 -openmp "
INC=$HOME/install/include
LIB=$HOME/install/lib

icc $FLAGS genData.c -o gen_$APP

icc $FLAGS ../train.c -I . -I $INC  -L $LIB -lnlopt -lm -o train_$APP.off

icc $FLAGS -Wno-unknown-pragmas -no-offload -O3 ../train.c -I . -I $INC \
				-L $LIB -lnlopt -lm -o train_$APP.omp

icc $FLAGS -Wno-unknown-pragmas -no-offload ../pred.c -I . -lm -o pred_$APP

FLAGS+=" -mmic -Wno-unknown-pragmas"
INC=$HOME/install_mic/include
LIB=$HOME/install_mic/lib

icc $FLAGS ../train.c -I . -I $INC   -L $LIB -lnlopt -lm -o train_$APP.mic

