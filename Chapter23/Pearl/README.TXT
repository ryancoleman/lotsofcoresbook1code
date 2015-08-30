Copyright (2012)2 (03-2014)3 Intel Corporation All Rights Reserved.

The source code contained or described herein and all documents related to the source code ("Material") are owned by Intel Corporation or its suppliers or licensors. Title to the Material remains with Intel Corporation or its suppliers and licensors. The Material may contain trade secrets and proprietary and confidential information of Intel Corporation and its suppliers and licensors, and is protected by worldwide copyright and trade secret laws and treaty provisions. No part of the Material may be used, copied, reproduced, modified, published, uploaded, posted, transmitted, distributed, or disclosed in any way without Intel’s prior express written permission. 
No license under any patent, copyright, trade secret or other intellectual property right is granted to or conferred upon you by disclosure or delivery of the Materials, either expressly, by implication, inducement, estoppel or otherwise. Any license under such intellectual property rights must be express and approved by Intel in writing.

Include any supplier copyright notices as supplier requires Intel to use.

Include supplier trademarks or logos as supplier requires Intel to use, preceded by an asterisk. An asterisked footnote can be added as follows: *Third Party trademarks are the property of their respective owners.

Unless otherwise agreed by Intel in writing, you may not remove or alter this notice or any other notice embedded in Materials by Intel or Intel’s suppliers or licensors in any way.”




iso-3dfd directory contains example of Intel(R) Xeon Phi(TM) implementation of

***************************************
16-th order 3D Finite Difference Stencil
***************************************

VERSION 1.0

Simple implementation for 3D FD stencil: 16-th order (49 points) with symmetric
and constant coefficients.

Questions: leonardo.borges@intel.com or philippe.thierry@intel.com or cedric.andreolli@intel.com


NOTE:
*************************
The results presented here are obtained with:
Intel Xeon Phi 7120 CO Coprocessor
2 x Intel E2697 V2 @ 2.7GHz


SETTINGS:
*************************
run_on_mic.pl:

Before trying to run the iso3DFD package on Xeon phi, you will need to setup your intel compiler path
in the file "run_on_mic.pl".

Xeon Phi modes:

On Xeon phi, it's possible to enable/disable ECC and to enable/disable turbo. You need to be root.
For turbo, you can use the command:
micsmc --turbo enable mic0
or
micsmc --turbo disable mic0

For ECC:
(ecc on)
micctrl -r
micctrl -w
micsmc --ecc enable mic0
micctrl -b
micctrl -w

(ecc off)
micctrl -r
micctrl -w
micsmc --ecc disable mic0
micctrl -b
micctrl -w


The fastest configuration for Xeon Phi is Ecc off and Turbo on.


NATIVE MODEL BUILD & RUN:
*************************


______________Quick Version________________

##If you want to run the latest version (dev09) on MIC:
./launch.sh


##If you want to run a specific version on the MIC:
./launch version=devXX


The available versions go from dev00 to dev09.

For example run:
./launch.sh version=dev08
n1=368 n2=390 n3=1300 nreps=100 num_threads=244 HALF_LENGTH=8
n1_thrd_block=368 n2_thrd_block=2 n3_thrd_block=26
allocating prev, next and vel: total 2135.19 Mbytes
-------------------------------
time:           2.89 sec
throughput:  5841.98 MPoints/s
flops:        356.36 GFlops

(Note that those results can be achieved with ECC off and Turbo On)




##For Xeon:
export BLOCK_LOOP="-DBLOCK_X_Z_Y";export KMP_AFFINITY="compact";export OPTI="-O3";export PREFETCH="-opt-prefetch-distance=88,22"

make clean version=dev07;make cleanall;make build version=dev07 model=cpu

./run_on_xeon.pl bin/iso3dfd_dev07_cpu_avx.exe 208 1227 1305 100 1216 49 74 compact 24



______________Detailled Version________________

##You can build a specific version of the code by using:
make build version=devXX

The available versions go from dev00 to dev09.


##You can also build a specific version for CPU or for MIC
make build version=devXX model=cpu

The default behavior is to build the program for CPU.


##You can also specify the vectorization level that you want to use:
make build version=devXX model=cpu simd=avx
