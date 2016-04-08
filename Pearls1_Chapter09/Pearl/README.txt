The initial version of this N-Body direct kernel is based on the description from the paper
"Test-driving Intel(r) Xeon Phi(tm) coprocessors with a basic N-body
simulation" by Andrey Vladimirov and Vadim Karpusenko 
(http://research.colfaxinternational.com/post/2013/01/07/Nbody-Xeon-Phi.aspx)

It relies on OpenMP and MKL. To compile and run just execute 'make' after
loading your compilation environment. By default it will try to execute on
mic0 but that can be adjusted by editing the Makefile.

You can run it also manually. There are two options:

* Using ssh (assumes that there is an NFS share between the host and the coprocessor):

	ssh mic0 LD_LIBRARY_PATH=$MIC_LD_LIBRARY_PATH $(pwd)/nbody

* Using micnativeloadex:

	SINK_LD_LIBRARY_PATH=$MIC_LD_LIBRARY_PATH /opt/intel/mic/bin/micnativeloadex ./nbody

Note that using micnativeloadex the last core of the coprocessor will not be
available to the application.

The application recognizes the following options:
	
	-s              Run only single precision version
	-d              Run only double precision version
	-i iters        Number of time steps to simulate
	-n particles    Number of particles to simulate

For example: 
	ssh mic0 LD_LIBRARY_PATH=$MIC_LD_LIBRARY_PATH $(pwd)/nbody -i 5 -n 20000 -s 
