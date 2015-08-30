# Note: gfortran > 4.4 is needed.  

  FC  := gfortran
  F90 := gfortran
  CC  := gcc

  ifdef MPI
  ifdef USE_MPI_WRAPPERS
    F90 = mpif90
  endif
  endif

  F90FLAGS += -J$(mdir) -I $(mdir)
  FFLAGS   += -J$(mdir) -I $(mdir)
  CFLAGS   += -std=c99 -Wall

  ifdef NDEBUG
    F90FLAGS += -O3
    FFLAGS   += -O3
    CFLAGS   += -O3
  else
    F90FLAGS += -g -O1 -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid -finit-real=nan
    FFLAGS   += -g -O1 -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid -finit-real=nan
    CFLAGS   += -g -O1
  endif

  ifdef OMP
    F90FLAGS += -fopenmp
    FFLAGS   += -fopenmp
    CFLAGS   += -fopenmp
  endif

  ifdef GPROF
    F90FLAGS += -pg
    FFLAGS   += -pg
    CFLAGS   += -pg
  endif
