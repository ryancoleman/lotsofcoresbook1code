# Makefile for MP-LABS (ZSC-3D-STD-MPI)
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
OBJS   = common.o main.o collision.o collisionrelax.o finaldump.o init.o   \
         memalloc.o parameters.o postcollision.o poststream.o relaxstats.o \
         stats.o stream.o vgrid.o vtkdump.o

#Link into an excutable
ZSC-3D-STD-MPI : $(OBJS); $(MPIFC) $(LFLAGS) $(OBJS) -o ZSC-3D-STD-MPI$(BIN_SUFFIX)

#Compile the common modules
common.o : common.F90; $(MPIFC) $(CFLAGS) common.F90

#Compile the subroutines
main.o           : main.F90 common.F90 collision.F90 collisionrelax.F90 finaldump.F90    \
                   init.F90 memalloc.F90 parameters.F90 postcollision.F90 poststream.F90 \
                   relaxstats.F90 stats.F90 stream.F90 vgrid.F90                       \
                   vtkdump.F90; $(MPIFC) $(CFLAGS) main.F90
collision.o      : common.F90 collision.F90; $(MPIFC) $(CFLAGS) collision.F90
collisionrelax.o : common.F90 collisionrelax.F90; $(MPIFC) $(CFLAGS) collisionrelax.F90
finaldump.o      : common.F90 finaldump.F90; $(MPIFC) $(CFLAGS) finaldump.F90
init.o           : common.F90 init.F90 stats.F90 vtkdump.F90; $(MPIFC) $(CFLAGS) init.F90
memalloc.o       : common.F90 memalloc.F90; $(MPIFC) $(CFLAGS) memalloc.F90
parameters.o     : common.F90 parameters.F90; $(MPIFC) $(CFLAGS) parameters.F90
postcollision.o  : common.F90 postcollision.F90; $(MPIFC) $(CFLAGS) postcollision.F90
poststream.o     : common.F90 poststream.F90; $(MPIFC) $(CFLAGS) poststream.F90
relaxstats.o     : common.F90 relaxstats.F90; $(MPIFC) $(CFLAGS) relaxstats.F90
stats.o          : common.F90 stats.F90; $(MPIFC) $(CFLAGS) stats.F90
stream.o         : common.F90 stream.F90; $(MPIFC) $(CFLAGS) stream.F90
vgrid.o          : common.F90 vgrid.F90; $(MPIFC) $(CFLAGS) vgrid.F90
vtkdump.o        : common.F90 vtkdump.F90; $(MPIFC) $(CFLAGS) vtkdump.F90

#Build options
all: ZSC-3D-STD-MPI

install:
	mkdir -p $(INSTALL_DIR)
	mv -v ./ZSC-3D-STD-MPI$(BIN_SUFFIX) $(INSTALL_DIR)

clean:
	rm -f ./*.o ./*.mod

build: ZSC-3D-STD-MPI

mrMonk:
	mkdir -p $(INSTALL_DIR)
	mv ./ZSC-3D-STD-MPI$(BIN_SUFFIX) $(INSTALL_DIR)
	rm -f ./*.o ./*.mod

