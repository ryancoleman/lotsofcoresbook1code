CXX=icpc

GIT_VERSION:=$(shell sh -c './GIT-VERSION-GEN')
COMPILER_VERSION:="$(CXX)-$(shell $(CXX) --version | head -n1 | cut -d' ' -f4)"

SSE_CXXFLAGS=-DSSE -xSSE4.2
KNC_CXXFLAGS=-DKNC -mmic
AVX_CXXFLAGS=-DAVX -xAVX
AVX2_CXXFLAGS=-DAVX -xCORE-AVX2

ARCH_CXXFLAGS=$(SSE_CXXFLAGS)

BASE_CXXFLAGS=-DGIT_VERSION=\"$(GIT_VERSION)\" -DCOMPILER_VERSION=\"$(COMPILER_VERSION)\" \
	      -DDOUBLE \
              -Wall -wd167 -ggdb \
	      -openmp \
              -O3 -ipo -restrict

all: run-tile compare

run-tile: run-tile.cpp pcl-hydro-core.cpp pcl-hydro-vcore.cpp pcl-hydro-util.cpp pcl-hydro-params.cpp vtkfile.cpp timeseries.cpp pcl-hydro.hpp arch.hpp Makefile
	$(CXX) -o $@ run-tile.cpp pcl-hydro-core.cpp pcl-hydro-vcore.cpp pcl-hydro-util.cpp pcl-hydro-params.cpp vtkfile.cpp timeseries.cpp $(BASE_CXXFLAGS) $(ARCH_CXXFLAGS)

compare: compare.cpp timeseries.cpp Makefile
	$(CXX) -o $@ compare.cpp timeseries.cpp $(BASE_CXXFLAGS)

clean:
	rm -rf run-tile compare
