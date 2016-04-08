########################################################################
#
# 'make' generates the executable 'printThreadCoreAssignment.x'
#
# How to run: e.g.
#
#    $> export MIC_ENV_PREFIX=MIC
#    $> export MIC_KMP_AFFINITY=compact
#    $> ./printThreadCoreAssignment.x
#
# NOTE: the default number of thread groups is 4, each of size 4.
#       that is, 16 threads are created on the Xeon Phi in the default case.
#       you can adjust these numbers by compiling via, e.g.,
#       'make GROUPS=12 GROUP_SIZE=7'
#
# Output: see sample 'output.log'
#
########################################################################
