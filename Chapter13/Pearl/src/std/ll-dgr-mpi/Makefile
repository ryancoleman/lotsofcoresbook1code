# Makefile for MP-LABS (LL-2D-DGR-MPI)
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
OBJS   = common.o finaldump.o hydrodynamics_g.o hydrodynamics_f.o init.o main.o\
         memalloc.o mpi_update_f.o mpi_update_g.o mpi_update_hydro.o           \
         mpi_update_rho_f.o mpi_update_rho_g.o parameters.o poststream_g.o     \
         poststream_f.o prestream_g.o prestream_f.o stats.o update_f.o         \
         update_g.o vgrid.o vtkplane.o

#Link into an excutable
LL-2D-DGR-MPI : $(OBJS); $(MPIFC) $(LFLAGS) $(OBJS) -o LL-2D-DGR-MPI$(BIN_SUFFIX)

#Compile the common modules
common.o : common.F90; $(MPIFC) $(CFLAGS) common.F90

#Compile the subroutines
main.o                 : common.F90 finaldump.F90 hydrodynamics_g.F90               \
                         hydrodynamics_f.F90 init.F90 memalloc.F90 mpi_update_f.F90   \
                         mpi_update_g.F90 mpi_update_hydro.F90 mpi_update_rho_f.F90 \
                         mpi_update_rho_g.F90 parameters.F90 poststream_g.F90       \
                         poststream_f.F90 prestream_g.F90 prestream_f.F90 stats.F90   \
                         update_f.F90 update_g.F90 vgrid.F90 vtkplane.F90             \
                         main.F90; $(MPIFC) $(CFLAGS) main.F90
finaldump.o            : common.F90 finaldump.F90; $(MPIFC) $(CFLAGS) finaldump.F90
hydrodynamics_g.o      : common.F90 hydrodynamics_g.F90; $(MPIFC) $(CFLAGS) hydrodynamics_g.F90
hydrodynamics_f.o      : common.F90 hydrodynamics_f.F90; $(MPIFC) $(CFLAGS) hydrodynamics_f.F90
init.o                 : common.F90 mpi_update_rho_g.F90 stats.F90 vtkplane.F90       \
                         init.F90; $(MPIFC) $(CFLAGS) init.F90
memalloc.o             : common.F90 memalloc.F90; $(MPIFC) $(CFLAGS) memalloc.F90
mpi_update_f.o         : common.F90 mpi_update_f.F90; $(MPIFC) $(CFLAGS) mpi_update_f.F90
mpi_update_g.o         : common.F90 mpi_update_g.F90; $(MPIFC) $(CFLAGS) mpi_update_g.F90
mpi_update_hydro.o     : common.F90 mpi_update_hydro.F90; $(MPIFC) $(CFLAGS) mpi_update_hydro.F90
mpi_update_rho_f.o     : common.F90 mpi_update_rho_f.F90; $(MPIFC) $(CFLAGS) mpi_update_rho_f.F90
mpi_update_rho_g.o     : common.F90 mpi_update_rho_g.F90; $(MPIFC) $(CFLAGS) mpi_update_rho_g.F90
parameters.o           : common.F90 parameters.F90; $(MPIFC) $(CFLAGS) parameters.F90
poststream_g.o         : common.F90 poststream_g.F90; $(MPIFC) $(CFLAGS) poststream_g.F90
poststream_f.o         : common.F90 poststream_f.F90; $(MPIFC) $(CFLAGS) poststream_f.F90
prestream_g.o          : common.F90 prestream_g.F90; $(MPIFC) $(CFLAGS) prestream_g.F90
prestream_f.o          : common.F90 prestream_f.F90; $(MPIFC) $(CFLAGS) prestream_f.F90
stats.o                : common.F90 stats.F90; $(MPIFC) $(CFLAGS) stats.F90
update_f.o             : common.F90 update_f.F90; $(MPIFC) $(CFLAGS) update_f.F90
update_g.o             : common.F90 update_g.F90; $(MPIFC) $(CFLAGS) update_g.F90
vgrid.o                : common.F90 vgrid.F90; $(MPIFC) $(CFLAGS) vgrid.F90
vtkplane.o             : common.F90 vtkplane.F90; $(MPIFC) $(CFLAGS) vtkplane.F90

#Build options
all: LL-2D-DGR-MPI

install:
	mkdir -p $(INSTALL_DIR)
	mv ./LL-2D-DGR-MPI$(BIN_SUFFIX) $(INSTALL_DIR)

clean:
	rm -f ./*.o ./*.mod

build: LL-2D-DGR-MPI

mrMonk:
	mkdir -p $(INSTALL_DIR)
	mv ./LL-2D-DGR-MPI$(BIN_SUFFIX) $(INSTALL_DIR)
	rm -f ./*.o ./*.mod

