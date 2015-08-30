ifdef TRACE
.PHONY: _trace _value
_trace: ; @$(MAKE) --no-print-directory TRACE= $(TRACE) ='$$(warning TRACE $(TRACE))$(shell $(MAKE) TRACE=$(TRACE)_value)'
_value: ; @echo '$(value $(TRACE))'
endif
OLD_SHELL := $(SHELL)
SHELL = $(warning Building $@$(if $<, (from $<))$(if $?, ($? newer)))$(OLD_SHELL) -x

include $(ROOT)/Makefile.include

.SUFFIXES:
.SUFFIXES: .f90 .F90 .o 

.PHONY: allobs
allobs:    $(OBJS)

%.o : %.f90
	@set -e; \cd $(<D); $(FC) $(FCFLAGS) -c $(<F)

%.o : %.F90
	@set -e; \cd $(<D); $(FC) $(FCCPP) $(FCFLAGS) -c $(<F)

