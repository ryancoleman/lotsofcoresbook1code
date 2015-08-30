# Makefile for MP-LABS (LBS3D-MPI)
#
# make CC=compiler CFLAGS="compiler flags" LFLAGS="linker flags" target
#
# Copyright 2013 Carlos Rosales Fernandez and The University of Texas at Austin.
# Copyright 2008 Carlos Rosales Fernandez, David S. Whyte and IHPC (A*STAR).
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

CFLAGS += $(CPP_FLAGS)
CFLAGS += $(OMP_FLAGS)
CFLAGS += $(ARCH_FLAGS)

LFLAGS += $(OMP_FLAGS)
LFLAGS += $(ARCH_FLAGS)

#Object files
OBJS = common.o main.o collision.o finalsave.o init.o memalloc.o parameters.o relaxstats.o postcollision.o poststream.o stats.o stream.o vgrid.o vtksave.o

#Link into an excutable
LBM : $(OBJS); $(MPIFC) $(LFLAGS) $(OBJS) -o lbs3d-mpi$(BIN_SUFFIX)

#Compile the common modules
common.o : common.F90; $(MPIFC) $(CFLAGS) common.F90

#Compile the subroutines
main.o           : main.F90 common.F90 collision.F90 finalsave.F90        \
                   init.F90 memalloc.F90 parameters.F90 postcollision.F90 \
                   poststream.F90 relaxstats.F90 stats.F90 stream.F90     \
                   vgrid.F90 vtksave.F90; $(MPIFC) $(CFLAGS) main.F90
collision.o      : common.F90 collision.F90;  $(MPIFC) $(CFLAGS) collision.F90
finalsave.o      : common.F90 finalsave.F90;  $(MPIFC) $(CFLAGS) finalsave.F90
init.o           : common.F90 init.F90;       $(MPIFC) $(CFLAGS) init.F90
memalloc.o       : common.F90 memalloc.F90;   $(MPIFC) $(CFLAGS) memalloc.F90
parameters.o     : common.F90 parameters.F90; $(MPIFC) $(CFLAGS) parameters.F90
relaxstats.o     : common.F90 relaxstats.F90; $(MPIFC) $(CFLAGS) relaxstats.F90
stats.o          : common.F90 stats.F90;      $(MPIFC) $(CFLAGS) stats.F90
stream.o         : common.F90 stream.F90;     $(MPIFC) $(CFLAGS) stream.F90
vgrid.o          : common.F90 vgrid.F90;      $(MPIFC) $(CFLAGS) vgrid.F90
vtksave.o        : common.F90 vtksave.F90;    $(MPIFC) $(CFLAGS) vtksave.F90
postcollision.o  : common.F90 postcollision.F90; $(MPIFC) $(CFLAGS) postcollision.F90
poststream.o     : common.F90 poststream.F90;    $(MPIFC) $(CFLAGS) poststream.F90

#Build options
build: LBM

install:
	mkdir -p $(INSTALL_DIR)
	mv -v ./lbs3d-mpi$(BIN_SUFFIX) $(INSTALL_DIR)

clean:
	-rm -f ./*.o ./*.mod

mrMonk:
	mkdir -p $(INSTALL_DIR)
	mv -v ./lbs3d-mpi$(BIN_SUFFIX) $(INSTALL_DIR)
	rm -f ./*.o ./*.mod



