OArchive file(s) in this directory dated in 2014 will correspond to 
version(s) used in the book.

The most current version can be retrieved from
https://github.com/WeiqunZhang/miniSMC

As noted in the book, "The best performance shown on Figure 4.19 for
the Intel Xeon Phi coprocessor is 12.99 seconds". Too late to modify
the manuscript before publication, the authors found that a few
additional Fortran compiler options increases performance to achieve
9.11 seconds for the coprocessor (1.43x improvement vs . Figure 4.19)
while keeping the processor performance the same. The compiler options
for the improved performance are:

-fno-alias -noprec-div -no-prec-sqrt -fimf-precision=low
 -fimf-domain-exclusion=15 -align array64byte -opt-assumesafe-padding
 -opt-streaming-stores always -opt-streaming-cache-evict=0.

These new best results on coprocessor use -1, 16, 16 for the thread
blocking specified in inputs_SMC file and SIMD directives are no
longer needed in kernels.F90 file when using -fno-alias compiler
option.  The chapter discusses Figure 4.19 as best results for
coprocessor.

The files here reflect all these changes for higher performance.

See the README in the Pearl directory to build and run.
