#
# $Id: GNUmakefile 24210 2013-05-11 19:49:27Z jhammond $
#

LIBRARY = libdriver.a

OBJ_OPTIMIZE = driver_input.o 

OBJ =       opt_drv.o oniom.o lbfgs.o lbfgs-b.o lbfgs-b_driver.o lbfgs_driver.o lbfgs_utils.o

USES_BLAS = opt_drv.F oniom.F lbfgs.F lbfgs-b.F lbfgs-b_driver.F

include ../config/makefile.h
include ../config/makelib.h

#FFLAGS  = -qsource -qxref=full -bloadmap:load.map 
