These are the sources and Makefile for the Morton ordering chapter.  Typing make on the command line will build 4 executables:

zmattr		- Morton order matrix transpose for Xeon
miczmattr	- Morton order matrix transpose for Xeon-Phi
zmatmul		- Morton order matrix multiply for Xeon
miczmatmul	- Morton order matrix multiply for Xeon-Phi



Execute zmattr and zmatmul from the command line and supply matrix dimension as the only command line parameter.  For example:

zmatmul 1024

will multiply two 1024x1024 matrices on Xeon.

Depending on your Xeon-Phi setup, copy miczmattr and miczmatmul to your working directory on Xeon-Phi and execute on the Xeon-Phi command line.

Support files:

zorder2d.c, zorder2d.h	- code to compute Morton order indices
timer.c, timer.h	- code used to time operations
