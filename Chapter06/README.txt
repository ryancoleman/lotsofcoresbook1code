Archive file(s) in this directory dated in 2014 will correspond to 
version(s) used in the book.

This version differs slightly from the source code in the book as it arranges
the computation in 'packets' of evaluations, which is the way the authors
constructed the production fault-tree evaluation code. As it turns out, this
helped with scaling on both the CPU and Xeon Phi.

Note: The performance of the chapter code on Xeon vs. Intel Xeon Phi is size
dependent as the Intel Xeon Phi needs to have enough work to keep the cores
and vector units busy. 

Following is one example where Intel Xeon Phi beats Intel Xeon

From the build directory on the host:

./fte  ../samples/sample3 loadedModel 1000000

From the home directory on a coprocessor:

./fte  ../samples/sample3 loadedModel 10000000




