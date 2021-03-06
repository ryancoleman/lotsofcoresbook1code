TARGET = pd.x
SRC = pd_main.cpp pd.cpp pd_force.cpp pd_update.cpp pd_energy.cpp
OBJ = $(SRC:%.cpp=%.o)

# particles are simulated within a 3D volume:
X=16 # x-extent
Y=16 # y-extent
Z=16 # z-extent

# number of thread groups for the computation
G=1

CXX = icpc
LD = icpc
CXXFLAGS = -O3 -xHost -fno-alias -opt-assume-safe-padding -openmp -DX=$(X) -DY=$(Y) -DZ=$(Z) -DNUM_GROUPS=$(G) -std=c++11 -Wall
LDFLAGS = -O3 -xHost -openmp -Wall

MATH_MODE = "fast2"

ifeq ($(MATH_MODE),"fast")
CXXFLAGS += -DFAST_MATH
endif

ifeq ($(MATH_MODE),"fast2")
CXXFLAGS += -DFAST_FAST_MATH
endif

all: $(TARGET)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

pd_force.o : pd_force.cpp
	$(CXX) $(CXXFLAGS) -vec-report=7 -c $< -o $@

pd_energy.o : pd_energy.cpp
	$(CXX) $(CXXFLAGS) -vec-report=7 -c $< -o $@

pd_update.o : pd_update.cpp
	$(CXX) $(CXXFLAGS) -vec-report=7 -c $< -o $@

$(TARGET) : $(OBJ)
	$(LD) $(LDFLAGS) -o $(TARGET) $(OBJ)

clean:
	rm -f *~ *.o *.optrpt $(TARGET)

run: $(TARGET)
	export output_X16_Y16_Z16_240threads_500PD_steps.log; export MIC_KMP_AFFINITY=compact,granularity=thread; ./$(TARGET)
