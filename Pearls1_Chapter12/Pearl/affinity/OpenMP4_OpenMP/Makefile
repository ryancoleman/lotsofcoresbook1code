TARGET = printThreadCoreAssignment.x
SRC = printThreadCoreAssignment.cpp
OBJ = $(SRC:%.cpp=%.o)

CXX = icpc
LD = icpc
CXXFLAGS = -O3 -openmp -std=c++11 -Wall
LDFLAGS = -O3 -openmp -Wall

GROUPS = 4
GROUP_SIZE = 4

all: $(TARGET)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -DGROUPS=$(GROUPS) -DGROUP_SIZE=$(GROUP_SIZE) -c $< -o $@

$(TARGET) : $(OBJ)
	$(LD) $(LDFLAGS) -o $(TARGET) $(OBJ)

clean:
	rm -f *.x *.o