
 QUICK START GUIDE
 Chapter 21:  High Performance Ray Tracing
 Building and Running The Tutorial On Intel Xeon Phi Coprocessors

 This tutorial is part of the Embree Kernel Library, which is included in the
 directory embree-source-2.3.3.  Alternatively, the latest version of Embree
 can be found at the following URL.  The tutorial is compiled as part of the
 Embree build process.

    https://github.com/embree

 What follows is a quick-start guide to building and running Embree on Linux
 hosts equipped with Intel Xeon Phi Coprocessors.  Please refer to the file
 README_Embree.txt for a full list of supported platforms, build dependencies,
 and complete build instructions.


 Linux Build Dependencies

   1. ISPC Compiler
      https://ispc.github.com
      export PATH=path-to-ispc:$PATH

   2. GLUT
      sudo yum install cmake.x86_64
      sudo yum install freeglut.x86_64 freeglut-devel.x86_64
      sudo yum install libXmu.x86_64 libXi.x86_64
      sudo yum install libXmu-devel.x86_64 libXi-devel.x86_64

   3. CMake


 Linux Build Instructions for Xeon Phi

   1. cd embree

   2. mkdir build

   3. cd build

   4. cmake ..

   5. ccmake ..

      The default build options should be correct in many cases, but verify
      that the XEON_ISA flag is set appropriately for your your system.  In
      addition, to build for Xeon Phi set the XEON_PHI_ISA flag to 'ON'.
      After doing so, press 'c' to configure, 'g' to generate the Make files,
      and 'q' to quit.

   6. cmake ..

   7. make


 Linux Execution Instructions for Xeon Phi

   1. cd embree/build

   2. export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH

      This path is used on the host to locate the Embree library and tutorial
      executable components built for Xeon.

   3. export SINK_LD_LIBRARY_PATH=$PWD:$SINK_LD_LIBRARY_PATH

      This path is used on the coprocessor to locate the Embree library and
      tutorial executable components built for Xeon Phi.

   4. ./tutorial00_xeonphi

      This tutorial runs on the host and offloads several compute kernels to
      the Xeon Phi coprocessor.  Please refer to README_Embree.txt for more
      information on command line and keyboard options that can be used with
      this tutorial, as well as a description of the other Embree tutorials.

