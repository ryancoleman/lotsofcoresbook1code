The codes here correspond to samples illustrating the "High Performance Parallelism Pearls" book chapter 11 "Dynamic load balancing using OpenMP 4.0".
The codes have been written by Gilles Civario and Michael Lysaght from the Irish Centre for High-End Computing and is released as is under MIT license. The text of the license should be available in this directory on the file license.txt.

The codes are a gradual evolution of an N-Body gravitational sample, from an usual OpenMP parallel version to a version able to take full advantage of all processors and co-processors available on the machine.
The various steps are as follow:
 - nbody-v0.cc: OpenMP as usual, able to run either on host Xeon processors or in native mode on one hosted Xeon Phi co-processor.
 - nbody-v1.cc: offload mode version of the code, running on the host and offloading the computations to one Xeon phi device.
 - nbody-v2.cc: uses all hosting Xeon cores and one Xeon Phi device at the same time. Load balancing between the two parts is dynamic and should always distribute the optimal amount of work to each sides of the machine. In this code, no overlapping is done during the data transfers between host and device.
 - nbody-v2.1.cc: same as previous with only an extra refinement to overlap sends and receives of data between host and device.
 - nbody-v3.cc: final version able of using an arbitrary number of devices on a machine, included none. This version will always use all Xeon and Xeon Phi co-processors it finds, and balance their respective workload as evenly as possible. Data transfers are reasonably overlapped but the code could be improved in that regards (if a performance difference between versions 2 and 2.1 showed it was beneficial).

The codes have been tested and are known to work with beta versions of the Intel C++ compiler 15.0. The GNU C++ compiler version 4.9.0 also permits to compile and run the codes, although it doesn't give access to any co-processor since the offloading directives are not honoured with this compiler version.
Earlier versions of either the Intel or the GNU compiler may or may not work, although later versions are likely to.

To compile the codes, adjust the compiler and compiler options in the Makefile. Samples for the Intel and GNU compilers are provided for reference.
Then just time "make".
All going well, this should generate binaries with the same base names as the corresponding source files, postfixed with "s" for single precision versions and "d" for double precision versions. If you chosen the Intel compiler, two extra native Xeon Phi binaries should have been compiled, with names postfixed with ".mic".

To run the various codes, just type "./nbody-vxp" with x being the version number and p being the precision. You can also adjust the number of threads to use by setting either OMP_NUM_THREADS or MIC_ENV_PREFIX to "MIC" and MIC_OMP_NUM_THREADS for the host and device respectively. That said, for the Intel compiler at least, not setting any OpenMP variables has proven the most effective since default values are just right.

For versions with offload mode enabled (all but version 0) when compiled with the Intel compiler, one might find the OFFLOAD_DEVICES environment variable useful for experimenting with different numbers or indexes of Xeon Phi co-processors, such as:
  > OFFLOAD_DEVICES="0,1" ./nbody-v3s     # run with devices 0 and 1 visible
  > OFFLOAD_DEVICES="1" ./nbody-v3s       # run with only the 2nd device visible
  > OFFLOAD_DEVICES="" ./nbody-v3s        # run with devices hidden

And finally, if compiled the Xeon Phi native versions of the initial code, you can run them the following way:
  > scp *.mic mic0:/tmp
  > ssh mic0 /tmp/*s.mic             # single precision version
  > ssh mic0 /tmp/*d.mic             # double precision version


