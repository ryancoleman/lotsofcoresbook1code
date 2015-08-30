README.txt
==========

In order to be able to use the OpenCL examples, you will need a working OpenCL implementation (usually available from your processor vendor's website) and a working C/C++ compiler that supports the latest C++11 standard. We've found that a gcc 4.8 or later works fine, as is a recent version of Intel's icc should also work well (v12 or later). You may also be able to get this working with clang.

To build the examples, make sure you have copied (and if necessary modified) one of the supplied "def" files to "make.def". For example, if working on Apple's OS X, copy apple.def to make.def. If working with an AMD OpenCL implementation on Linux, then copy linux_amd.def to make.def, etc. Example def files are supplied for: on Linux: Intel, AMD and Nvidia; on Windows: Intel; and for Mac OS X.

If you system has multiple OpenCL platforms and devices, make sure you select your desired target device type in your make.def *before* building the OpenCL host program. The three most likely options are:

To target a CPU (your OpenCL platform must support this - e.g. Intel or AMD):

  DEVICE = CL_DEVICE_TYPE_CPU

To target a GPU (your OpenCL platform must support this - e.g. Intel, AMD or Nvidia):

  DEVICE = CL_DEVICE_TYPE_GPU

To target an accelerator, such as an Intel Xeon Phi coprocessor:

  DEVICE = CL_DEVICE_TYPE_ACCELERATOR

If you're not sure what devices are supported by your OpenCL implementation, you can try:

  DEVICE = CL_DEVICE_TYPE_DEFAULT

Once you have an appropriate make.def for your system and OpenCL installation,  you should be able to type "make" (nmake if using Windows) to build the OpenCL host executable, in this case "matmul".

You can then run the executable from the command line. On Linux this would be with a simple command such as "./matmul"

