#$Id: GNUmakefile 19708 2010-10-29 18:04:21Z d3y133 $

include ../config/makefile.h

        OBJ = \
	      mp2_grad.o \
              mp2_energy.o \
	      mp2_grad_doit.o \
	      mp2_make_cphf_rhs.o \
	      mp2_make_tot_dens.o \
	      mp2_make_tot_wdens.o \
	      mp2_make_dip_mom.o \
              mp2_pstat.o \
	      mp2_input.o \
              mp2_make_piy.o mp2_make_nos.o

OBJ_OPTIMIZE = mp2_back_transform.o \
               mp2_read_tiajb.o \
               mp2_pijab_uhf.o \
               mp2_wijab_uhf.o \
               mp2_lai_uhf.o \
               mp2_lai_fock_uhf.o mp2_memory.o

ifneq (,$(findstring SGI, $(TARGET)))
   OBJ += mp2_make_tuhf.o
else
   OBJ_OPTIMIZE += mp2_make_tuhf.o
endif
ifeq ($(TARGET),DECOSF)
	FC = f77 -notransform_loops
#  -O5 alone breaks on mp2_make_tuhf.o ; -notransform loops keeps O5 going
# EA 4/20/200
endif


USES_BLAS = mp2_energy.F mp2_lai_uhf.F \
            mp2_make_tuhf.F mp2_pijab_uhf.F \
            mp2_wai_uhf.F mp2_wijab_uhf.F mp2_back_transform.F \
            mp2_read_tiajb.F

    LIBRARY = libmp2.a


LIB_DEFINES =  -DDEBUG_PRINT


# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */


include ../config/makelib.h


