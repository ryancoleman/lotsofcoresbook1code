# Authors: Gilles Civario and Michael Lysaght
# Copyright 2014 Irish Centre for High-End Computing under MIT license
# See license.txt for the full license of this software

# Intel compiler (tested with version 15.0.0)
CXX=icpc
CXXFLAGS=-O3 -fopenmp -fp-model fast=2
HOST=-xhost
DEVICE=-mmic -openmp-link=static

# or GNU compiler (tested with version 4.9.0)
#CXX=g++
#CXXFLAGS=-O3 -fopenmp -ffast-math
#HOST=-march=native -mtune=native
#DEVICE=

SRC:=$(wildcard *.cc)
EXE:=$(patsubst %.cc,%s,$(SRC)) $(patsubst %.cc,%d,$(SRC))
ifeq "$(CXX)" "icpc"
   EXEMIC=nbody-v0s.mic nbody-v0d.mic
endif

all: $(EXE) $(EXEMIC)

%s: %.cc
	$(CXX) $(CXXFLAGS) $(HOST) -o $@ $<

%d: %.cc
	$(CXX) $(CXXFLAGS) $(HOST) -DDOUBLE -o $@ $<

%s.mic: %.cc
	$(CXX) $(CXXFLAGS) $(DEVICE) -o $@ $<

%d.mic: %.cc
	$(CXX) $(CXXFLAGS) $(DEVICE) -DDOUBLE -o $@ $<

clean:
	rm -f $(EXE) $(EXEMIC)


