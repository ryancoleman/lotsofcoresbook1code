    F90 := ifort
    FC  := ifort
    CC  := icc

    ifdef MPI
    ifdef USE_MPI_WRAPPERS
      F90 = mpiifort
    endif
    endif

    FFLAGS   = -module $(mdir) -I $(mdir) -fno-alias  -align array64byte
    F90FLAGS = -module $(mdir) -I $(mdir) -fno-alias  -align array64byte
    CFLAGS   = -std=c99

    ifdef OMP
      FFLAGS   += -openmp -openmp-report2
      F90FLAGS += -openmp -openmp-report2
      CFLAGS   += -openmp -openmp-report2
    endif

    ifdef MIC
      FFLAGS   += -mmic -no-prec-div -no-prec-sqrt -fimf-precision=low -fimf-domain-exclusion=15 -qopt-assume-safe-padding -opt-streaming-stores always -opt-streaming-cache-evict=0
      F90FLAGS += -mmic -no-prec-div -no-prec-sqrt -fimf-precision=low -fimf-domain-exclusion=15 -qopt-assume-safe-padding -opt-streaming-stores always -opt-streaming-cache-evict=0
      CFLAGS   += -mmic
    endif

    ifdef NDEBUG
      F90FLAGS += -O3
      FFLAGS   += -O3
      CFLAGS   += -O3
    else
      F90FLAGS += -g -traceback -O0 #-check all -warn all -u 
      FFLAGS   += -g -traceback -O0 #-check all -warn all -u 
      CFLAGS   += -g -Wcheck
    endif

    ifdef GPROF
      FFLAGS   += -p
      F90FLAGS += -p
      CFLAGS   += -p
    endif
