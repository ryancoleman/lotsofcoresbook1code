##beacon
MPILIB= -L$(MPICH_HOME)/lib -lmpi
MKLDIR=$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
BLASLIB=-lmkl_blas95_lp64 -lpthread -lm
#BLASLIB=-lmkl_blas95_lp64 -lmkl_blacs_intelmpi_lp64 -lpthread -lm
LAPACKLIB=-lmkl_lapack95_lp64
PBLAS_DIR=
#PBLAS_DIR=-I$(MKLROOT)/include
#SCALAPACKLIB= -lmkl_scalapack_lp64 -lmkl_intel_lp64 
SCALAPACKLIB=../ooc_solver/scalapack/libscalapack.a 

PORTLIB = pdlaprnt.f pzlaprnt.f pslaprnt.f pclaprnt.f

##beacon
LIBS = \
    -L$(MKLDIR) $(LAPACKLIB) $(BLASLIB) $(SCALAPACKLIB) 

# CFLAGS_O = -O2 -profile-functions
CFLAGS_O = -O2 
CFLAGS_g = -g
CFLAGS = -DADD_ -us $(CFLAGS_g) \
    $(PBLAS_DIR)  \
    -DUSE_MIC=1 \
    -mkl -offload-option,mic,compiler,"-z defs" \
    -DUSE_PROFILE=1 \
    -DUSE_PROF=1 \
    -Wno-deprecated

##beacon
CC=mpiicc

# --------------------
# options for gfortran
# --------------------
FFLAGS_g = -g
FFLAGS_O =  -O2 
FFLAGS = $(FFLAGS_O) \
    -DUSE_MIC \
    -mkl -offload-option,mic,compiler,"-z defs" \
    -DUSE_PROF \
    -DUSE_PROFILE

##beacon
FC = mpiifort
