TARGET = dataTransfer.x
SRC = dataTransfer.cpp
OBJ = $(SRC:%.cpp=%.o)

CXX = mpiicpc
LD = mpiicpc
CXXFLAGS = -O3 -openmp -std=c++11 -Wall
LDFLAGS = -O3 -openmp -lrt -Wall

all: $(TARGET)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET) : $(OBJ)
	$(LD) $(LDFLAGS) -o $(TARGET) $(OBJ)

clean:
	rm -f *.x *.o