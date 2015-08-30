# $Id: GNUmakefile 19696 2010-10-29 16:53:42Z d3y133 $


               OBJ = sp1111.o
      OBJ_OPTIMIZE = genr70.o init70.o qout70.o \
                     sinfo.o sp0000.o sp0001.o \
                     sp0011.o sp0101.o sp0111.o \
                     tabgen.o 

    LIBRARY = libnwints.a


  USES_BLAS = sp1111.F

include ../../config/makefile.h
include ../../config/makelib.h

$(OBJ):	auxvar.fh geom.fh maxc.fh  pgeom.fh qgeom.fh \
	const.fh  ginf.fh miscg.fh picon.fh shllfo.fh


