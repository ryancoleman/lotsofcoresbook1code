# Minimum GCC4.7 is required to be installed on the machine for compilation with Intel(r) ComposerXE
# GCC4.9 is supported only for compilation for Intel(r) Xeon(r)

ifeq ($(COMP),gcc)
  OPENMP_FLAGS:=-openmp
else
  COMP:=icc
  OPENMP_FLAGS:=-openmp -par-affinity="granularity=fine,scatter"
endif

ifeq ($(PLAT),knc)
  ARCH:= mic
  TBB_LIB:=$(TBBROOT)/lib/mic/gcc4.4
  NUMA_LIB:=
  ifneq ($(TEAMS),)
    TEAM_DEF:=-DMIC_TEAMS=$(TEAMS)
  else
    TEAM_DEF:=-DMIC_TEAMS=20
  endif
  ifeq ($(TASKS),)
    TASK_DEF:=-D_NUM_TASKS=20
  else
    TASK_DEF:=-D_NUM_TASKS=$(TASKS)
  endif
  NUMA_CMD:=micnativeloadex ./gems
  FLAT_CMD:=micnativeloadex ./gems -a flat
else
  ARCH:= avx
  NUMA_LIB:=-lnuma
  TBB_LIB:= $(TBBROOT)/lib/intel64/gcc4.4
  NUMA_CMD:=./gems
  FLAT_CMD:=./gems flat
  ifeq ($(TASKS),)
    TASK_DEF:=-D_NUM_TASKS=60
  else
    TASK_DEF:=-D_NUM_TASKS=$(TASKS)
  endif
endif

COPT+=-O3 -m$(ARCH) -g
COPT+=-D__ENABLE_TRAPPING__ -D__ENABLE_AFFINITY__

RUN_DEPENDS:=run-flat
ifeq ($(BASELINE),TBB)
  COPT+=-D__TBB_BASELINE__
  RUN_DEPENDS:=run-numa
endif

ifeq ($(BASELINE),OPENMP)
  COPT+=-D__OMP_BASELINE__ $(OPENMP_FLAGS)
  RUN_DEPENDS:=run-numa
endif

ifeq ($(BASELINE),CILK)
  COPT+=-D__CILK_BASELINE__
  RUN_DEPENDS:=run-numa
endif

ifeq ($(BASELINE),TASK)
  COPT+=-D__TASK_BASELINE__
  RUN_DEPENDS:=run-numa
endif

run: gems $(RUN_DEPENDS)

run-flat: run-numa
	$(FLAT_CMD)

run-numa:
	$(NUMA_CMD)

gems: Gems_Code.cpp arenas.h observers.h ittnotify_static.o Makefile build.cmd
	$(COMP) -o gems -std=c++11 $(COPT) -I$(TBBROOT)/include Gems_Code.cpp ittnotify_static.o -L$(TBB_LIB) -ltbb -ltbbmalloc -lpthread $(NUMA_LIB) ${TEAM_DEF} $(TASK_DEF)

ittnotify_static.o: Makefile
	$(COMP) -c -m$(ARCH) -O2 ../ittnotify/ittnotify_static.c -o ittnotify_static.o

clean:
	rm gems *.o

