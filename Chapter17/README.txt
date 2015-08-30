Archive file(s) in this directory dated in 2014 will correspond to
version(s) used in the book, which was downloaded from
http://www.nwchem-sw.org/download.php?f=Nwchem-6.5.revision26243-src.2014-09-10.tar.gz

The source code can be downloaded from http://www.nwchem-sw.org/download.php.

This version was downloaded from wget
http://www.nwchem-sw.org/images/Nwchem-6.5.xeonphi.tar.gz

Building:
   cd nwchem-6.5;
   sh build4xeonphi.sh


Running:
   cd examples
   mpirun -np 4 ../bin/LINUX64/nwchem uracil_dimer.nw

The variable NWC_RANKS_PER_DEVICE can be used to turn off the Xeon Phi
offloading.  If you set NWC_RANKS_PER_DEVICE=0 then NWChem will not
run on the Xeon Phi, just on the Xeon CPUs.

the examples directory contains an input file (uracil_dimer.nw and a
few output files).

The output files are from a 16 cores Sandy Bridge system equipped with
two KNC cards (if you are interested in output from runs using
multiples nodes in a cluster, let me know).  The only env. variable
you might want to set is OMP_NUM_THREADS (that will be used for the
Xeon CPU threading), since threading, affinity, etc .. for Xeon Phi is
automatically set by the code (but can be still turned off with
NWC_RANKS_PER_DEVICE as shown above). Two processes are automatically
assigned for offloading to each KNC card, therefore the remaining
process will work on the Xeon CPUs instead.

----- Additional Information about NWChem ---------


Compilation instruction at
http://www.nwchem-sw.org/index.php/Compiling_NWChem

Documentation at
http://www.nwchem-sw.org/index.php/Release65:NWChem_Documentation

Release notes at
http://www.nwchem-sw.org/index.php/NWChem_6.5


