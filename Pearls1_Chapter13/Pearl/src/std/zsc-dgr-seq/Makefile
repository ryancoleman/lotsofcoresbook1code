# Makefile for MP-LABS (ZSC-2D-DGR)
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
OBJS   = common.o collision_f.o collision_g.o differentials.o finaldump.o \
         init.o main.o memalloc.o parameters.o stats.o stream.o update.o  \
         vtkplane.o

#Link into an excutable
ZSC-2D-DGR : $(OBJS); $(FC) $(LFLAGS) $(OBJS) -o ZSC-2D-DGR$(BIN_SUFFIX)

#Compile the common modules
common.o : common.F90; $(FC) $(CFLAGS) common.F90

#Compile the subroutines
main.o          : common.F90 collision_f.F90 collision_g.F90 differentials.F90        \
                  finaldump.F90 init.F90 memalloc.F90 parameters.F90 stats.F90 stream.F90 \
                  update.F90 vtkplane.F90 main.F90 ; $(FC) $(CFLAGS) main.F90
collision_f.o   : common.F90 collision_f.F90; $(FC) $(CFLAGS) collision_f.F90
collision_g.o   : common.F90 collision_g.F90; $(FC) $(CFLAGS) collision_g.F90
differentials.o : common.F90 differentials.F90; $(FC) $(CFLAGS) differentials.F90
finaldump.o     : common.F90 finaldump.F90; $(FC) $(CFLAGS) finaldump.F90
init.o          : common.F90 differentials.F90 stats.F90 vtkplane.F90                 \
                  init.F90; $(FC) $(CFLAGS) init.F90
memalloc.o      : common.F90 memalloc.F90; $(FC) $(CFLAGS) memalloc.F90
parameters.o    : common.F90 parameters.F90; $(FC) $(CFLAGS) parameters.F90
stats.o         : common.F90 stats.F90; $(FC) $(CFLAGS) stats.F90
stream.o        : common.F90 stream.F90; $(FC) $(CFLAGS) stream.F90
update.o        : common.F90 differentials.F90; $(FC) $(CFLAGS) update.F90
vtkplane.o      : common.F90 vtkplane.F90; $(FC) $(CFLAGS) vtkplane.F90

# Make options
all: ZSC-2D-DGR

install:
	mkdir -p $(INSTALL_DIR)
	mv -v ./ZSC-2D-DGR$(BIN_SUFFIX) $(INSTALL_DIR)

clean:
	rm -f ./*.o ./*.mod

build: ZSC-2D-DGR

mrMonk:
	mkdir -p $(INSTALL_DIR)
	mv ./ZSC-2D-DGR$(BIN_SUFFIX) $(INSTALL_DIR)
	rm -f ./*.o ./*.mod
