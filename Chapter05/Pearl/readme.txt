There are two make files in this folder

Makefile             intended to be run on Linux with "make"
MakefileWindows      intended to be run on Windows with "nmake"

To build both Xeon Phi and host programs on Linux

   make all

To build both Xeon Phi and host programs on Windows

   nmake -f MakefileWindows all

When building just host programs:

    change the target all to taget all_host

When building just Xeon Phi programs

    change target from all to all_xphi

There are additional make targets, please consult the appropriate make file.

There are two Linux bash shell scripts and one Windows batch file

run_host.sh, run_xphi.sh and run_host.bat

The run_xphi.sh is intended to be run on the MIC under its Linux O/S

The run_host.sh runs all the host programs and writes the statistics to the terminal window.
*** prior to running issue

   export MyNumberOfCores=(your number of cores here)

The above specifies the number of cores on your system.

The scripts require HT be enabled (as this is a demonstration of HT enabled programs).

On Windows

   SET NUMBER_OF_CORES=(your number of cores here)

The run_xphi.sh script is configured for a 60 core 240 thread coprocessor.
You may have to edit this script should your coprocessor differ.

When running the bash shell scripts, include the -x option such that
you obtain the name of the program being run in addition to its output.

  bash -x ./run_xphi.sh

The above issued in an ssh session to the MIC

  bash -x ./run_host.sh

The above issued on the host Linux terminal window.

  run_host

The above issued on the Windows CMD window.

All builds and runs on host are assumed to be initiated in the src folder
(wherever you place that)

*** It is a requirement that the proper environment be set up for your Intel C++ x64 compiler
and when Xeon Phi are being built, that a current version of MPSS be installed.

Jim Dempsey