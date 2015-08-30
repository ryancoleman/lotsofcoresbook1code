TARGET = pd.x
SRC = pd_main.cpp pd.cpp pd_force.cpp pd_update.cpp pd_energy.cpp
OBJ = $(SRC:%.cpp=%.o)

# particles are simulated within a 3D volume:

X=16 # x-extent
Y=16 # y-extent
Z=16 # z-extent

# number of threads for the computation
T=32

CXX = icpc
LD = icpc
CXXFLAGS = -O3 -xHost -fno-alias -opt-assume-safe-padding -openmp -DX=$(X) -DY=$(Y) -DZ=$(Z) -DNUM_THREADS=$(T) -std=c++11 -Wall
LDFLAGS = -O3 -xHost -openmp -Wall

MATH_MODE = "fast2"

ifeq ($(MATH_MODE),"fast")
CXXFLAGS += -DFAST_MATH
endif

ifeq ($(MATH_MODE),"fast2")
CXXFLAGS += -DFAST_FAST_MATH
endif

PINNING_SCHEME=scatter
ifeq ($(T),240)
PINNING_SCHEME=compact
endif

all: $(TARGET)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

pd_force.o : pd_force.cpp
	$(CXX) $(CXXFLAGS) -vec-report=7 -c $< -o $@

pd_update.o : pd_update.cpp
	$(CXX) $(CXXFLAGS) -vec-report=7 -c $< -o $@

pd_energy.o : pd_energy.cpp
	$(CXX) $(CXXFLAGS) -vec-report=7 -c $< -o $@

$(TARGET) : $(OBJ)
	$(LD) $(LDFLAGS) -o $(TARGET) $(OBJ)

clean:
	rm -f *.o $(TARGET)

run: $(TARGET)
	export MIC_ENV_PREFIX=MIC; export MIC_KMP_AFFINITY=$(PINNING_SCHEME),granularity=thread; ./$(TARGET)
