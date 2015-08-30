# Makefile for MP-LABS (LL-2D-STD)
#
# make FC=compiler CFLAGS="compiler flags" LFLAGS="linker flags" target
#
# Copyright 2008 Carlos Rosales Fernandez, David S. Whyte and IHPC (A*STAR).
# Copyright 2013 Carlos Rosales Fernandez and The University of Texas at Austin.
#
# This file is part of MP-LABS.
#
# MP-LABS is free software: you can redistribute it and/or modify it under the
# terms of the GNU GPL version 3 or (at your option) any later version.
#
# MP-LABS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# MP-LABS, in the file COPYING.txt. If not, see <http://www.gnu.org/licenses/>.


SHELL = /bin/sh

#Object files
OBJS   = common.o finaldump.o hydrodynamics.o init.o main.o memalloc.o \
         parameters.o poststream.o prestream.o stats.o vtkplane.o

#Link into an excutable
LL-2D-STD : $(OBJS); $(FC) $(LFLAGS) $(OBJS) -o LL-2D-STD$(BIN_SUFFIX)

#Compile the common modules
common.o : common.F90; $(FC) $(CFLAGS) common.F90

#Compile the subroutines
main.o               : common.F90 finaldump.F90 hydrodynamics.F90 init.F90 memalloc.F90 \
                       parameters.F90 poststream.F90 prestream.F90 stats.F90          \
                       vtkplane.F90 main.F90; $(FC) $(CFLAGS) main.F90
finaldump.o          : common.F90 finaldump.F90; $(FC) $(CFLAGS) finaldump.F90
hydrodynamics.o      : common.F90 hydrodynamics.F90; $(FC) $(CFLAGS) hydrodynamics.F90
init.o               : common.F90 stats.F90 vtkplane.F90 init.F90; $(FC) $(CFLAGS) init.F90
memalloc.o           : common.F90 memalloc.F90; $(FC) $(CFLAGS) memalloc.F90
parameters.o         : common.F90 parameters.F90; $(FC) $(CFLAGS) parameters.F90
poststream.o         : common.F90 poststream.F90; $(FC) $(CFLAGS) poststream.F90
prestream.o          : common.F90 prestream.F90; $(FC) $(CFLAGS) prestream.F90
stats.o              : common.F90 stats.F90; $(FC) $(CFLAGS) stats.F90
vtkplane.o           : common.F90 vtkplane.F90; $(FC) $(CFLAGS) vtkplane.F90


#Build options
all: LL-2D-STD

install:
	mkdir -p $(INSTALL_DIR)
	mv -v ./LL-2D-STD$(BIN_SUFFIX) $(INSTALL_DIR)

clean:
	rm -f ./*.o ./*.mod

build: LL-2D-STD

mrMonk:
	mkdir -p $(INSTALL_DIR)
	mv ./LL-2D-STD$(BIN_SUFFIX) ../../../bin
	rm -f ./*.o ./*.mod

