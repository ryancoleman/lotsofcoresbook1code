
Code Sample for Data Transfer Using the Intel COI Library

* Benchmark Performance Disclaimer

The source code in this directory are example source code used in the book chapter Data Transfer Using the 
Intel COI Library. The examples are meant to demotrate how different COI buffer types can be used and their
relative performance. Depending on the particular applications and system configurations, there are further 
optimizations can be done to improve the performance. On some systems the result might be worse. The benchmark 
results from these examples should not be used for any purpose other than learning the Intel COI library.

* How to Build and Run

This code was tested with GCC 4.4.7 and ICC 14.0 on RHEL 6.3.

Ensure the Intel MPSS stack 3.2 or later are installed on your system:
https://software.intel.com/en-us/articles/intel-manycore-platform-software-stack-mpss

Also make sure the COI SDK is installed as part of the MPSS installation, check in the following directories:

/opt/intel/mic/coi
/usr/include/intel-coi/

If you want to use Intel ICC compiler, please edit the Makefile to make the nessasory changes (see comments
in the Makefile).

Build:
- make

Run:
- cd release
- ./source_stream
- ./source_normal
- ./source_pinned

Note that each source program will iterate three times, where the first run will warm up the system.