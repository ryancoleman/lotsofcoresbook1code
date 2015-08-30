This is a direct N-body code that runs on Xeon and Xeon Phi.
To run on Xeon go to the directory "avx" and type "make" to compile and run the code.

To run on Xeon Phi go to the directory "mic".
First set the variables "LIB" and "MIC" in the Makefile.
Then typing "make" will compile and run the code.

"make build" will only compile.
"make run" will only run.
"make clean" will remove the object file and executable.
"make" defaults to "make all", which will run "make clean", "make build", then "make run".
