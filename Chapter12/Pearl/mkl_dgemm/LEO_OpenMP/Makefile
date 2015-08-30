TARGET = mkl_dgemm.x
SRC = mkl_dgemm.cpp
OBJ = $(SRC:%.cpp=%.o)

CXX = icpc
LD = icpc
CXXFLAGS = -O3 -openmp -mkl -std=c++11 -Wall
LDFLAGS = -O3 -openmp -mkl -lrt -Wall

all: $(TARGET)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET) : $(OBJ)
	$(LD) $(LDFLAGS) -o $(TARGET) $(OBJ)

clean:
	rm -f *.x *.o