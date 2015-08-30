# Makefile for MP-LABS (LBS3D-DEVEL-OFFLOAD)
#
# make FC=compiler CFLAGS="compiler flags" LFLAGS="linker flags" target
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
CFLAGS += -offload-option,mic,compiler,"-O3"

LFLAGS += $(OMP_FLAGS)
LFLAGS += $(ARCH_FLAGS)
LFLAGS += -offload-option,mic,ld,"-O3"

#Object files
OBJS   = common.o main.o collision.o collisionmic.o finalsave.o init.o initmic.o   \
         memalloc.o parameters.o relaxstats.o postcollision.o postcollisionmic.o   \
         poststream.o poststreammic.o stats.o stream.o streammic.o packcpu.o       \
         packcpu_f.o packmic.o packmic_f.o updatephi.o updatephimic.o vtksave.o

#Link into an excutable
LBM : $(OBJS); $(FC) $(LFLAGS) $(OBJS) -o lbs3d-devel-offload$(BIN_SUFFIX)

#Compile the common modules
common.o : common.F90; $(FC) $(CFLAGS) common.F90

#Compile the subroutines
main.o           : main.F90 common.F90 collision.F90 collisionmic.F90 finalsave.F90     \
                   init.F90 initmic.F90 memalloc.F90 parameters.F90 postcollision.F90   \
                   postcollisionmic.F90 poststream.F90 poststreammic.F90 relaxstats.F90 \
                   stats.F90 stream.F90 streammic.F90 packcpu.F90 packmic.F90           \
                   packcpu_f.F90 packmic_f.F90 updatephi.F90 updatephimic.F90           \
                   vtksave.F90; $(FC) $(CFLAGS) main.F90
collision.o        : common.F90 collision.F90;     $(FC) $(CFLAGS) collision.F90
collisionmic.o     : common.F90 collisionmic.F90;  $(FC) $(CFLAGS) collisionmic.F90
finalsave.o        : common.F90 finalsave.F90;     $(FC) $(CFLAGS) finalsave.F90
init.o             : common.F90 init.F90;          $(FC) $(CFLAGS) init.F90
initmic.o          : common.F90 initmic.F90;       $(FC) $(CFLAGS) initmic.F90
memalloc.o         : common.F90 memalloc.F90;      $(FC) $(CFLAGS) memalloc.F90
parameters.o       : common.F90 parameters.F90;    $(FC) $(CFLAGS) parameters.F90
relaxstats.o       : common.F90 relaxstats.F90;    $(FC) $(CFLAGS) relaxstats.F90
stats.o            : common.F90 stats.F90;         $(FC) $(CFLAGS) stats.F90
stream.o           : common.F90 stream.F90;        $(FC) $(CFLAGS) stream.F90
streammic.o        : common.F90 streammic.F90;     $(FC) $(CFLAGS) streammic.F90
vtksave.o          : common.F90 vtksave.F90;       $(FC) $(CFLAGS) vtksave.F90
packcpu.o          : common.F90 packcpu.F90;       $(FC) $(CFLAGS) packcpu.F90
packmic.o          : common.F90 packmic.F90;       $(FC) $(CFLAGS) packmic.F90
packcpu_f.o        : common.F90 packcpu_f.F90;     $(FC) $(CFLAGS) packcpu_f.F90
packmic_f.o        : common.F90 packmic_f.F90;     $(FC) $(CFLAGS) packmic_f.F90
updatephi.o        : common.F90 updatephi.F90;     $(FC) $(CFLAGS) updatephi.F90
updatephimic.o     : common.F90 updatephimic.F90;  $(FC) $(CFLAGS) updatephimic.F90
postcollision.o    : common.F90 postcollision.F90;    $(FC) $(CFLAGS) postcollision.F90
postcollisionmic.o : common.F90 postcollisionmic.F90; $(FC) $(CFLAGS) postcollisionmic.F90
poststream.o       : common.F90 poststream.F90;       $(FC) $(CFLAGS) poststream.F90
poststreammic.o    : common.F90 poststreammic.F90;    $(FC) $(CFLAGS) poststreammic.F90

#Build options
build: LBM

install:
	mkdir -p $(INSTALL_DIR)
	mv -v ./lbs3d-devel-offload$(BIN_SUFFIX) $(INSTALL_DIR)

clean:
	-rm -f ./*.o ./*.mod

mrMonk:
	mkdir -p $(INSTALL_DIR)
	mv -v ./lbs3d-devel-offload$(BIN_SUFFIX) $(INSTALL_DIR)
	rm -f ./*.o ./*.mod

