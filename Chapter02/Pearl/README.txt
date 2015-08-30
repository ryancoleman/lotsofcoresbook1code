pcl-hydro, a reimplementation of the github hydro2d code.

contact: jason.sewall@intel.com

LICENSE
---------

See LICENSE.txt

COMPILING
---------

To compile for a given architecture, edit the Makefile line ARCH_CXXFLAGS= to be set to one of the 4 variables above it. You can also leave it empty for a (data-)serial version.

RUNNING
-------

To run it, do

$ ./run-tile -i test.nml

which will run based on the parameters in the test.nml file included with the program.

You can override keyword parameters on the command line; I recommend overriding nx and ny to increase the spatial size of the problem on faster systems.

i.e.
$ ./run-tile -i test.nml --nx=1024 --ny=1024

If you need it to run for longer or shorter, you can change the --nstepmax and/or --tend parameters

The KNC version is native, so copy input files and use micnativeloadx or whatever you like.

THREADING
---------

Threading is determined by OpenMP environment variables.

OUTPUT
------

The program will print out a line for each step with some information about the timelevel/dt. These are good canaries in the coal mine; dt should be in the ~1e-3 range.

You can quiet this down by giving the --quiet or -q options; multiples of these will make it extra quiet.

The program will print its own configuration and timing info to stderr upon completion.

VALIDATION
----------

Giving it the --output <dirs>/<prefix> option will cause all <dirs>/ to be created (if needed) and the files <dirs>/<prefix>.idx and <dirs>/<prefix><seq>.pak to be created with one or more frames of output. Repeated invocations will cause truncation of existing files.

You can control the frequency of the output with --output-interval.

You can compare with another output with the compare program:

$ ./compare <path-to-idx1> <path-to-idx2>

which will print out norms between each frames of the datasets. The datasets must share dimensions, and you can restrict the frame ranges as you like.

MISC
----

The tiling might fail, especially for KNC, for small problems. This is because I insist on having at least SIMD_WIDTH cells in each tile in X and don't allow for other cases. This isn't a fundemental problem, but one I didn't take the time to address.
