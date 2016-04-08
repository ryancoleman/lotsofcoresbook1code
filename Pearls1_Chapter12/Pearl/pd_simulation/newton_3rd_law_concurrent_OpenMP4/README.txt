########################################################################
#
# 'make' generates the executable 'pd.x'
#
# How to run: e.g.
#
#    $> ./pd.x
#
# or use the Makefile (you can specifiy the number of
# thread groups and volume extent):
#
#    $> make run
#    $> make clean
#    $> make G=8 X=16 Y=32 Z=16 run
#
# Code was tested on Xeon Phi 3xxx, 5xxx, 7xxx with:
#
#    mpss 2.1, mpss 3.1, mpss 3.4
#    intel compilers 14.0.3, 15.0.0
#
########################################################################
