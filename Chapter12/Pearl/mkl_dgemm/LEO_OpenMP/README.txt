########################################################################
#
# 'make' generates the executable 'mkl_dgemm.x'
#
# How to run: e.g.
#
#    $> ./mkl_dgemm.x
#
# or use the bash script 'run.sh' to run our benchmarks
# (it takes long to complete)
#
#    $> ./run.sh
#
# Output: see sample output in 
#    mkl_dgemm_performance_size2048L_groupSize16_copyFraction0.0_xeonPhi7xxx_mpss3.1_compilerintel14.0.3_schedSetAffinity_xxx.log
#
# Code was tested on Xeon Phi 3xxx, 5xxx, 7xxx with:
#
#    mpss 2.1, mpss 3.1, mpss 3.4
#    intel compilers 14.0.3, 15.0.0
#
########################################################################
