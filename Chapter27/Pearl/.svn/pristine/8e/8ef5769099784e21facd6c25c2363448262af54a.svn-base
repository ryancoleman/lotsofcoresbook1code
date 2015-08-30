
#some advanced dependencies require user input on library locations
#(for e.g. the CSB interface (scheme 16) and the MKL interface (scheme 18).
#
#search for `#start' in this makefile on where to supply this info.

.SILENT:help

SLFLAGS=
DLFLAGS=
SCFLAGS=${CFLAGS} -static
DCFLAGS=${CFLAGS} -fPIC
LIBNAME=sparselib
CURRENTVERSION=1.6.0RC1
SHORTCURRENTVERSION=${CURRENTVERSION:%.=%}

FILE_SRC=FileToVT.hpp FileToVT.cpp mmio.cpp mmio.h
UTIL_SRC=util.cpp util.hpp MachineInfo.cpp MachineInfo.hpp
HILBERT_SRC=Matrix2HilbertCoordinates.cpp Matrix2HilbertCoordinates.hpp
BLOCKORDER_SRC=BlockOrderer.hpp SBDTree.cpp SBDTree.hpp MinCCS.hpp

FILE_C=$(filter %.cpp,${FILE_SRC})
UTIL_C=$(filter %.cpp,${UTIL_SRC})
HILBERT_C=$(filter %.cpp,${HILBERT_SRC})
BLOCKORDER_C=$(filter %.cpp,${BLOCKORDER_SRC})

FILE_DEPS=${FILE_SRC} ${FILE_C:%.cpp=%.o}
UTIL_DEPS=${UTIL_SRC} ${UTIL_C:%.cpp=%.o}
HILBERT_DEPS=${HILBERT_SRC} ${HILBERT_C:%.cpp=%.o}
BLOCKORDER_DEPS=${BLOCKORDER_SRC} ${BLOCKORDER_C:%.cpp=%.o}
DEPS=${FILE_DEPS} ${UTIL_DEPS} ${HILBERT_DEPS} ${BLOCKORDER_DEPS}

CILK=cilk/bin/cilk++

LFLAGS=-lpthread -lrt -lnuma

#start compiler choices
#1: GNU Compiler Collection (GCC)
CC=g++
# start different CFLAGS options
# a: debug mode
#CFLAGS=-ansi -Wall -pedantic -O0 -g -fprefetch-loop-arrays -funroll-loops -fopenmp
#CMSG="Debug flags"
# b: profiling mode
#CFLAGS=-ansi -Wall -pedantic -g -fprefetch-loop-arrays -funroll-loops -fopenmp -DNDEBUG -pg
#CMSG="Profiling flags"
# c: performance mode
CFLAGS=-ansi -Wall -pedantic -g -O3 -fprefetch-loop-arrays -funroll-loops -ffast-math -fno-exceptions -fopenmp -DNDEBUG
CMSG="Performance flags"
# end different CFLAGS options
#2: Intel Compiler Collection (ICC)
#IPATH=/opt/intel/composer_xe_2013.1.117/
# a: performance mode
#CC=icc -diag-disable 654 -diag-disable 68
#CILK=icc
#CFLAGS=-xHost -O3 -openmp -DNDEBUG
#CMSG="ICC performance flags"
# b: performance mode, Xeon Phi target
#CC=icc -diag-disable 654 -diag-disable 68
#CILK=icc
#CFLAGS=-mmic -O3 -openmp -DNDEBUG -D_NO_LIBNUMA
#CMSG="ICC performance flags -- MIC target"
#LFLAGS=-pthread -lrt
#end compiler choices

#start MKL section
#MKL_PATH=${IPATH}mkl/
#to use MKL, check if the above path is correct and un-comment the below
#MKL_INCLUDE_PATH=${MKL_PATH}include/
#CFLAGS:=${CFLAGS} -DWITH_MKL -I${MKL_INCLUDE_PATH}
#uncomment for Xeon targets:
#MKL_LIBRARY_PATH=${MKL_PATH}lib/intel64/
#uncomment for Xeon Phi target:
#MKL_LIBRARY_PATH=${MKL_PATH}lib/mic/
#un-comment first & second line in case of GCC, third & fourth line in case of ICC
#SLFLAGS:=${SLFLAGS} -Wl,--start-group ${MKL_LIBRARY_PATH}libmkl_intel_lp64.a ${MKL_LIBRARY_PATH}libmkl_gnu_thread.a ${MKL_LIBRARY_PATH}libmkl_core.a -Wl,--end-group
#DLFLAGS:=${DLFLAGS} -L${MKL_LIBRARY_PATH} -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core
#SLFLAGS:=${SLFLAGS} -Wl,--start-group ${MKL_LIBRARY_PATH}libmkl_intel_lp64.a ${MKL_LIBRARY_PATH}libmkl_intel_thread.a ${MKL_LIBRARY_PATH}libmkl_core.a -Wl,--end-group -ldl -liomp5
#DLFLAGS:=${DLFLAGS} -L${MKL_LIBRARY_PATH} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
#end of MKL section

#start CUDA section; to use CuSparse, uncomment the below lines
#CUDA_PATH=/usr/local/cuda-5.5/
#CFLAGS:=${CFLAGS} -DWITH_CUDA -I${CUDA_PATH}include/
#DLFLAGS:=${DLFLAGS} -L${CUDA_PATH}lib64/ -lcusparse -lcudart
#SLFLAGS:=${SLFLAGS} -L${CUDA_PATH}lib64/ -lcusparse -lcudart #CUDA only offers dynamic linkage for its libraries
#DEPS:=${DEPS} CuHyb.o
#end CUDA section

CLEAN=$(filter %.o,${DEPS})
CLEAN:=${CLEAN:%.o=%.static.o} ${CLEAN}

all: msg mm2htrp driver-all sbd2trp sbdmv mtx2crstrp McShared McDMV
	echo "Done!"

msg:
	echo "*** Compiling using ${CMSG} ***"

help:
	echo "Sparse library, version ${CURRENTVERSION}"
	echo ""
	echo "DEFAULT TARGETS (make all)"
	echo "=========================="
	echo "The default target makes, next to the library (see next section), the following utilities:"
	echo "   1. 'driver'        Parallel and sequential SpMV benchmark utility, statically linked to the sparse library"
	echo "   2. 'driver-shared' Same as 1, but using shared linkage"
	echo "   3. 'sbdmv'         Sequential SpMV benchmark according to the cache-oblivious (CO) 2D SBD strategy"
	echo "   4. 'McDMV'         Parallel explicit 2D SpMV benchmark augmented with CO 2D SBD (see above)"
	echo "   5. 'McShared'      Parallel SpMV benchmark according to the implicit synchronising 2D approach"
	echo "   6. 'mm2htrp'       Converts a  .mtx file into a binary Hilbert-ordered .trp file"
	echo "   7. 'sbd2trp'       Converts an .emm file, containing a matrix in 2D SBD format, into an SBD-ordered .trp file"
	echo "   8. 'mtx2crstrp'    Converts a  .mtx file to a CRS-ordered .trp file"
	echo ""
	echo "Execute each utility without arguments for details."
	echo ""
	echo "Note that some schemes make use of the MachineInfo utilility to determine the number of threads to run on;"
	echo "to manually define this number P, store P in the file \"hardware.info\" in the current directory, e.g.:"
	echo "    echo \"4\" > hardware.info; ./driver <program args>"
	echo "will do an SpMV benchmark using 4 threads if driver calls a parallel scheme."
	echo ""
	echo "LIBRARY COMPILATION (make all/libs)"
	echo "==================================="
	echo "'make libs'   attempts to compile both a static and shared library."
	echo "'make shared' Only compiles the shared library (.so)"
	echo "'make static' Only compiles the static library (.a)"
	echo ""
	echo "('make all' implies 'make libs')"
	echo "The libraries (.a / .so) are put into the ./lib/ directory."
	echo "The header files required for the library are put into ./include/"
	echo ""
	echo "OTHER MAKE TARGETS:"
	echo "==================="
	echo "'make McDirect'      - Experimental parallel explicit 2D SpMV multiplication (based on McShared, superceded by McDMV)."
	echo "'make intel'         - Does 'make all libs' using the intel compiler (icc) targeting intel architectures."
	echo "'make core2'         - Does 'make all libs' using the intel compiler (icc) targeting core2 architectures."
	echo "'make docs'          - Compiles javadoc documentation; will be generated in the ./doxy_out/ directory."
	echo ""
	echo "CLEANUP TARGETS:"
	echo "================"
	echo "'make clean'     - Cleans all .o files."
	echo "'make veryclean' - Implies 'make clean', also removes all compiled executables, the library directories './lib' and './include', and the doxygen directory."
	echo ""

intel:
	make CC=icc CFLAGS="-O2 -DNDEBUG -fprefetch-loop-arrays -funroll-loops -ansi -fPIC -m64 -msse3 -lrt" all libs

core2:
	make CC=icc CFLAGS="-O2 -DNDEBUG -fprefetch-loop-arrays -funroll-loops -ansi -fPIC -m64 -march=core2 -mtune=core2 -lrt" all libs

McDirect: McDirect.o static
	${CC} ${CFLAGS} -o McDirect McDirect.o lib/lib${LIBNAME}${CURRENTVERSION}.a -lpthread -lrt

McShared.o : McShared.cpp
	${CC} ${CFLAGS} -DTRIPLET_META -DTRIPLET_META_TYPE='unsigned int' -c -o McShared.o McShared.cpp

McShared.static.o : McShared.cpp
	${CC} ${SCFLAGS} -DTRIPLET_META -DTRIPLET_META_TYPE='unsigned int' -c -o McShared.static.o McShared.cpp

McShared: McShared.o static
	echo "...building McShared";
	${CC} ${CFLAGS} -DTRIPLET_META -DTRIPLET_META_TYPE='unsigned int' -o McShared McShared.o lib/lib${LIBNAME}${CURRENTVERSION}.a -lpthread -lrt

McDMV: McDMV.cpp static
	echo "...building McDMV";
	${CC} ${CFLAGS} -DTRIPLET_META -DTRIPLET_META_TYPE='unsigned int' -o $@ ${^:static=} lib/lib${LIBNAME}${CURRENTVERSION}.a ${LFLAGS}

libs: clean
	make CFLAGS="${CFLAGS} -DNDEBUG" static shared

cilk-static: static
	echo "Forcing WITH_CSB flag"
	echo "Forcing CILKPARALLEL flag"
	${CILK} -DWITH_CSB -DCILKPARALLEL ${CFLAGS} -o cilk-static driver.cpp ${FILE_C} ${UTIL_C} ${HILBERT_C} ${LFLAGS}

driver: static driver.static.o
	echo "...building driver application (statically linked)";
	${CC} ${CFLAGS} -fopenmp -o driver driver.static.o lib/lib${LIBNAME}${CURRENTVERSION}.a ${LFLAGS} ${SLFLAGS}

driver-shared: shared driver.o
	echo "...building driver application (dynamically linked)";
	${CC} ${CFLAGS} -fopenmp -o driver-shared driver.o -Llib/ lib/lib${LIBNAME}.so.${CURRENTVERSION} ${LFLAGS} ${DLFLAGS}

driver-all:
	make driver
	make driver-shared

micbench: static micbench.static.o
	${CC} ${CFLAGS} -o micbench micbench.static.o lib/lib${LIBNAME}${CURRENTVERSION}.a ${LFLAGS} ${SLFLAGS}

%.static.o: %.c %.h
	${CC} ${SCFLAGS} -c -o $@ ${^:%.h=}

%.static.o: %.cpp %.h
	${CC} ${SCFLAGS} -c -o $@ ${^:%.h=}

%.static.o: %.cpp %.hpp
	${CC} ${SCFLAGS} -c -o $@ ${^:%.hpp=}

%.static.o: %.cpp
	${CC} ${SCFLAGS} -c -o $@ ${^:%.h=}

%.o: %.c %.h
	${CC} ${DCFLAGS} -c -o $@ ${^:%.h=}

%.o: %.cpp %.h
	${CC} ${DCFLAGS} -c -o $@ ${^:%.h=}

%.o: %.cpp %.hpp
	${CC} ${DCFLAGS} -c -o $@ ${^:%.hpp=}

%.o: %.cpp
	${CC} ${DCFLAGS} -c -o $@ $^

static: ${DEPS:%.o=%.static.o}
	echo "...building static library";
	mkdir -p include || true
	cp *.hpp include/
	cp *.h include/
	mkdir -p lib || true;
	ar rcs lib${LIBNAME}${CURRENTVERSION}.a $(filter %.o,$^)
	mv lib${LIBNAME}${CURRENTVERSION}.a ./lib/

shared: ${DEPS}
	echo "...building dynamic library";
	mkdir -p include || true
	cp *.hpp include/
	cp *.h include/
	mkdir -p lib || true;
	${CC} 	${CFLAGS} \
		-shared -Wl,-soname,lib${LIBNAME}.so.${SHORTCURRENTVERSION} \
		-o ./lib/lib${LIBNAME}.so.${CURRENTVERSION} $(filter %.o,$^) -lc

docs: doxy.conf
	doxygen doxy.conf

mm2htrp: mm2htrp.cpp
	echo "...building mm2htrp";
	${CC} -Wall -O3 -ansi -o mm2htrp mm2htrp.cpp Matrix2HilbertCoordinates.cpp FileToVT.cpp mmio.cpp -lrt

sbd2trp: sbd2trp-main.cpp ${FILE_DEPS} Triplet.hpp ${BLOCKORDER_DEPS} ${UTIL_DEPS}
	echo "...building sbd2trp";
	${CC} ${CFLAGS} -o sbd2trp sbd2trp-main.cpp $(filter %.o,$^)

sbdmv: sbdmv.cpp ${FILE_DEPS} Triplet.hpp ${HILBERT_DEPS} ${UTIL_DEPS} ${BLOCKORDER_DEPS}
	echo "...building sbdmv";
	${CC} ${CFLAGS} -o sbdmv sbdmv.cpp $(filter %.o,$^) -lrt

mtx2crstrp: mtx2crstrp.cpp ${FILE_DEPS}
	echo "...building mtx2crstrp";
	${CC} ${CFLAGS} -o $@ mtx2crstrp.cpp $(filter %.o,$^)

clean:
	rm -f driver*.o
	rm -f micbench*.o
	rm -f McShared.o
	rm -rf $(filter %.o,${CLEAN})

veryclean: clean
	rm -f mm2htrp
	rm -f sbd2trp
	rm -f McDirect
	rm -f McShared
	rm -f McDMV
	rm -f mtx2crstrp
	rm -rf sbdmv
	rm -rf include
	rm -rf lib
	rm -rf doxy_out
	rm -f driver
	rm -f driver-shared

