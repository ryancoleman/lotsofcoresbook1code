#!/usr/bin/perl
# Copyright (2012)2 (03-2014)3 Intel Corporation All Rights Reserved. 
# The source code contained or described herein and all documents related 
# to the source code ("Material") are owned by Intel Corporation or its suppliers 
# or licensors. Title to the Material remains with Intel Corporation or its 
# suppliers and licensors. The Material contains trade secrets and proprietary 
# and confidential information of Intel or its suppliers and licensors. The 
# Material is protected by worldwide copyright and trade secret laws and treaty 
# provisions. No part of the Material may be used, copied, reproduced, modified, 
# published, uploaded, posted, transmitted, distributed, or disclosed in any way 
# without Intel’s prior express written permission.
#
# No license under any patent, copyright, trade secret or other intellectual property 
# right is granted to or conferred upon you by disclosure or delivery of the 
# Materials, either expressly, by implication, inducement, estoppel or otherwise. 
# Any license under such intellectual property rights must be express 
# and approved by Intel in writing.

# Script created by cedric.andreolli@intel.com

use strict;
use warnings;

#Setup your compiler path here !
my $compilerPath = "/opt/intel/composer_xe_2013_sp1.2.144";





my $argNumber = 12;

if(@ARGV < $argNumber){
	print "This script must be called with the following parameters:\n";
	print "\t-executable_name: The mic executable name\n";
	print "\t-micIP: The mic ip (mic0, mic1, etc.)\n";
	print "\t-n1: N1\n";
	print "\t-n2: N2\n"; 
	print "\t-n3: N3\n"; 
	print "\t-nb_iter: The number of iterations\n"; 
	print "\t-n1_block: The size of the n1_block\n"; 
	print "\t-n2_block: The size of the n2 block\n"; 
	print "\t-n3_block: The size of the n3 block\n"; 
	print "\t-kmp_affinity: The thread partitionning\n"; 
	print "\t-nb_cores: The number of cores\n"; 
	print "\t-nb_threads_by_core: The number of threads by core\n"; 
	exit();
}

my $exe = $ARGV[0];

my $micIP = $ARGV[1];
my $n1 = $ARGV[2];
my $n2 = $ARGV[3];
my $n3 = $ARGV[4];
my $nbIter = $ARGV[5];
my $n1Block = $ARGV[6];
my $n2Block = $ARGV[7];
my $n3Block = $ARGV[8];
my $kmpAffinity = $ARGV[9];
my $nbCores = $ARGV[10];
my $nbThreadsByCore = $ARGV[11];


my $nbThreads = $nbCores * $nbThreadsByCore;

print `scp $compilerPath/compiler/lib/mic/libiomp5.so $micIP:/tmp`;
print `scp $exe $micIP:/tmp`;

if($exe =~/\/.*\.exe/){
	($exe) = $exe =~/\/([^\/]*\.exe)/
}


print `ssh $micIP "cd /tmp; export LD_LIBRARY_PATH=.;ulimit -s 8192;export KMP_AFFINITY=$kmpAffinity;./$exe $n1 $n2 $n3 $nbThreads $nbIter $n1Block $n2Block $n3Block"`;
sleep(2);
print `ssh $micIP "rm -rf /tmp/$exe /tmp/libiomp5.so"`;
