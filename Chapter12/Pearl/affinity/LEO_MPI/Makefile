TARGET = printThreadCoreAssignment.x
SRC = printThreadCoreAssignment.cpp
OBJ = $(SRC:%.cpp=%.o)

GROUP_SIZE=4

CXX = mpiicpc
LD = mpiicpc
CXXFLAGS = -O3 -openmp -std=c++11 -DGROUP_SIZE=$(GROUP_SIZE) -Wall
LDFLAGS = -O3 -openmp -Wall


all: $(TARGET)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET) : $(OBJ)
	$(LD) $(LDFLAGS) -o $(TARGET) $(OBJ)

clean:
	rm -f *.x *.o