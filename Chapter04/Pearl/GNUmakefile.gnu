#MPI := 
MPI := t

#OMP := 
OMP := t

NDEBUG := t
#NDEBUG := 

### Edit compiler options in comps/$(COMP).mak
COMP := GNU
# COMP := Intel

K_USE_AUTOMATIC := t

### 
ifdef MPI
  # Set USE_MPI_WRAPPERS := t to use mpif90 or mpiifort, 
  # otherwise you need to specify mpi_include_dir, mpi_lib_dir, and mpi_libraries.
  USE_MPI_WRAPPERS := t
  ifndef USE_MPI_WRAPPERS
    MPI_HOME := /usr/local
    mpi_include_dir = $(strip $(MPI_HOME))/include
    mpi_lib_dir = $(strip $(MPI_HOME))/lib
    mpi_libraries += -lmpich -lpthread
  endif
endif

MKVERBOSE := t

#########################################################################

COMP := $(strip $(COMP))

ifdef MPI
  mpi_suffix 	:= .mpi
endif
ifdef OMP
  omp_suffix 	:= .omp
endif
ifndef NDEBUG
  debug_suffix 	:= .debug
endif

suf=$(COMP)$(debug_suffix)$(mpi_suffix)$(omp_suffix)

odir=.
mdir=.
tdir=.

tdir = t/$(suf)
odir = $(tdir)/o
mdir = $(tdir)/m

F_C_LINK := UNDERSCORE

SOURCE_HOME := src

fsources    := $(notdir $(wildcard $(SOURCE_HOME)/*.f))
f90sources  := $(notdir $(wildcard $(SOURCE_HOME)/*.f90))
F90sources  := $(notdir $(wildcard $(SOURCE_HOME)/*.F90))
csources    := $(notdir $(wildcard $(SOURCE_HOME)/*.c))

ifdef MPI
  f90sources := $(filter-out parallel_stubs.f90,$(f90sources))
else
  f90sources := $(filter-out parallel.f90,$(f90sources))
endif

ifdef OMP
  f90sources := $(filter-out omp_stubs.f90,$(f90sources))
else
  f90sources := $(filter-out omp.f90,$(f90sources))
endif

libraries =

fpp_flags =
fld_flags =

ifdef K_USE_AUTOMATIC
  F90PPFLAGS += -DK_USE_AUTOMATIC
else
  F90PPFLAGS =
endif

ifeq ($(wildcard comps/$(COMP).mak),)
   $(error "comps/$(COMP).mak does not exist")   
else 
  include comps/$(COMP).mak
endif

ifdef mpi_include_dir
  fpp_flags += -I $(mpi_include_dir)
endif

ifdef mpi_lib_dir
  fld_flags += -L $(mpi_lib_dir)
endif

MODDEP  :=  scripts/moddep.pl
MKDEP   :=  scripts/mkdep.pl

FPPFLAGS += $(fpp_flags) -I $(SOURCE_HOME)
LDFLAGS  += $(fld_flags) 
libraries += $(mpi_libraries)

CPPFLAGS += -DBL_FORT_USE_$(F_C_LINK) -I $(SOURCE_HOME)

objects = $(addprefix $(odir)/,       \
	$(sort $(f90sources:.f90=.o)) \
	$(sort $(F90sources:.F90=.o)) \
	$(sort $(fsources:.f=.o))     \
	$(sort $(csources:.c=.o))     \
	)
sources =                     \
	$(sort $(f90sources)) \
	$(sort $(F90sources)) \
	$(sort $(fsources)  ) \
	$(sort $(csources)  )

c_includes = $(SOURCE_HOME)
f_includes = $(SOURCE_HOME)

vpath %.c $(SOURCE_HOME)
vpath %.f $(SOURCE_HOME)
vpath %.f90 $(SOURCE_HOME)
vpath %.F90 $(SOURCE_HOME)

COMPILE.f   = $(FC)  $(FFLAGS) $(FPPFLAGS) -c
COMPILE.f90 = $(F90) $(F90FLAGS) $(FPPFLAGS) -c
COMPILE.F90 = $(F90) $(F90FLAGS) $(F90PPFLAGS) $(FPPFLAGS) -c

LINK.f90    = $(F90) $(F90FLAGS) $(FPPFLAGS) $(LDFLAGS)

default: main.$(suf).exe

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

clean::
	$(RM) ./*.o ./*.mod $(mdir)/*.mod $(odir)/*.o *.$(suf).exe *~
	$(RM) $(tdir)/f90.depends $(tdir)/c.depends
	$(RM) TAGS tags

realclean:: clean
	$(RM) -fr t
	$(RM) *.exe

TAGS:	$(sources)
	ctags -e --verbose=yes --fortran-kinds=+i $(abspath $^)

tags:	$(sources)
	ctags --verbose=yes --fortran-kinds=+i $^

# should prevent deletion of .o files
.SECONDARY: $(objects)


%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o main.$(suf).exe $< $(objects) $(libraries)
else
	@echo "Linking $@ ..."
	@$(LINK.f90) -o main.$(suf).exe $< $(objects) $(libraries)
endif

${odir}/%.o: %.f
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
ifdef MKVERBOSE
	$(COMPILE.f) -o $@ $<
else
	@echo "Building $< ..."
	@$(COMPILE.f) -o $@ $<
endif

${odir}/%.o: %.f90
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
ifdef MKVERBOSE
	$(COMPILE.f90) -o $@ $<
else
	@echo "Building $< ..."
	@$(COMPILE.f90) -o $@ $<
endif

${odir}/%.o: %.F90
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
ifdef MKVERBOSE
	$(COMPILE.F90) -o $@ $<
else
	@echo "Building $< ..."
	@$(COMPILE.F90) -o $@ $<
endif

${odir}/%.o: %.c
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
ifdef MKVERBOSE
	$(COMPILE.c) -o $@ $<
else
	@echo "Building $< ..."
	@$(COMPILE.c) -o $@ $<
endif

$(tdir)/f90.depends: $(fsources) $(f90sources) $(F90sources)
	@if [ ! -d $(tdir) ]; then mkdir -p $(tdir); fi
ifdef MKVERBOSE
	perl $(MODDEP) $(f_includes) --odir $(odir)  $^ > $(tdir)/f90.depends 
else
	@echo "Building f90/f dependency File ..."
	@perl $(MODDEP) $(f_includes) --odir $(odir) $^ > $(tdir)/f90.depends 
endif

$(tdir)/c.depends:  $(csources)
	@if [ ! -d $(tdir) ]; then mkdir -p $(tdir); fi
ifdef MKVERBOSE
	perl $(MKDEP) $(c_includes) --odir $(odir) $^ > $(tdir)/c.depends 
else
	@echo "Building c dependency File ..."
	@perl $(MKDEP) $(c_includes) --odir $(odir) $^ > $(tdir)/c.depends 
endif

ifneq ($(MAKECMDGOALS),realclean)
ifneq ($(MAKECMDGOALS),clean)
include $(tdir)/f90.depends

ifdef csources
include $(tdir)/c.depends
endif
endif
endif

#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)

#-----------------------------------------------------------------------------
# If a file requires special compiler flags (say -O0), one can do
# something like below, where SPECIALF90FLAGS is modified based on the
# result of "make print-F90FLAGS".
#
#SPECIALF90FLAGS := -O0 -Jt/GNU.mpi.omp/m -I t/GNU.mpi.omp/m -ftree-vectorize -fno-range-check -fopenmp
#${odir}/xxx.o: src/xxx.f90
#	$(F90) $(SPECIALF90FLAGS) $(FPPFLAGS) -c -o $@ $<
