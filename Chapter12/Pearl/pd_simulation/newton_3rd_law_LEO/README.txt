########################################################################
#
# 'make' generates the executable 'pd.x'
#
# How to run: e.g.
#
#    $> ./pd.x
#
# or use the Makefile (you can specifiy thread count and volume extent):
#
#    $> make run
#    $> make clean
#    $> make T=120 X=16 Y=32 Z=16 run
#
# Output: see sample output
#    output_X16_Y16_Z16_240threads_500PD_steps.log
#
# NOTE: for performance measurements comment '#define MEASUREMENT' 
# in 'pd_main.cpp'
#
# Code was tested on Xeon Phi 3xxx, 5xxx, 7xxx with:
#
#    mpss 2.1, mpss 3.1, mpss 3.4
#    intel compilers 14.0.3, 15.0.0
#
########################################################################
