# Makefile for MP-LABS (LL-2D-STD-MPI)
#
# make MPIFC=compiler CFLAGS="compiler flags" LFLAGS="linker flags" target
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
         parameters.o poststream.o prestream.o stats.o mpi_updatefg.o  \
        mpi_updaterho.o vgrid.o vtkplane.o

#Link into an excutable
LL-2D-STD-MPI : $(OBJS); $(MPIFC) $(LFLAGS) $(OBJS) -o LL-2D-STD-MPI$(BIN_SUFFIX)

#Compile the common modules
common.o : common.F90; $(MPIFC) $(CFLAGS) common.F90

#Compile the subroutines
main.o               : common.F90 finaldump.F90 hydrodynamics.F90 init.F90 memalloc.F90 \
                       parameters.F90 poststream.F90 prestream.F90 stats.F90          \
                       mpi_updatefg.F90 mpi_updaterho.F90 vgrid.F90 vtkplane.F90      \
                       main.F90; $(MPIFC) $(CFLAGS) main.F90
finaldump.o          : common.F90 finaldump.F90; $(MPIFC) $(CFLAGS) finaldump.F90
hydrodynamics.o      : common.F90 hydrodynamics.F90; $(MPIFC) $(CFLAGS) hydrodynamics.F90
init.o               : common.F90 stats.F90 vtkplane.F90 init.F90; $(MPIFC) $(CFLAGS) init.F90
memalloc.o           : common.F90 memalloc.F90; $(MPIFC) $(CFLAGS) memalloc.F90
parameters.o         : common.F90 parameters.F90; $(MPIFC) $(CFLAGS) parameters.F90
poststream.o         : common.F90 poststream.F90; $(MPIFC) $(CFLAGS) poststream.F90
prestream.o          : common.F90 prestream.F90; $(MPIFC) $(CFLAGS) prestream.F90
stats.o              : common.F90 stats.F90; $(MPIFC) $(CFLAGS) stats.F90
mpi_updatefg.o       : common.F90 mpi_updatefg.F90; $(MPIFC) $(CFLAGS) mpi_updatefg.F90
mpi_updaterho.o      : common.F90 mpi_updaterho.F90; $(MPIFC) $(CFLAGS) mpi_updaterho.F90
vgrid.o              : common.F90 vgrid.F90; $(MPIFC) $(CFLAGS) vgrid.F90
vtkplane.o           : common.F90 vtkplane.F90; $(MPIFC) $(CFLAGS) vtkplane.F90

#Build oprions
all: LL-2D-STD-MPI

install:
	mkdir -p $(INSTALL_DIR)
	mv -v ./LL-2D-STD-MPI$(BIN_SUFFIX) $(INSTALL_DIR)

clean:
	rm -f ./*.o ./*.mod

build: LL-2D-STD-MPI

mrMonk:
	mkdir -p $(INSTALL_DIR)
	mv ./LL-2D-STD-MPI$(BIN_SUFFIX) $(INSTALL_DIR)
	rm -f ./*.o ./*.mod
