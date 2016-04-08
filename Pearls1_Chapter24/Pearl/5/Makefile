CXX = icpc
CXXFLAGS = -g -O3 -fopenmp -DMETHOD=5
MICFLAGS  = 
HOSTFLAGS = 

CPUOBJECTS = \
	Transpose-CPU.o\
	Main-CPU.o

MICOBJECTS = \
	Transpose-MIC.o\
	Main-MIC.o

TARGETS = runme-CPU runme-MIC

%-CPU.o: %.cc
	$(CXX) $(CXXFLAGS) $(HOSTFLAGS) -DDOUBLE -c -o $@ $^

%-MIC.o: %.cc
	$(CXX) $(CXXFLAGS) $(MICFLAGS) -mmic -DDOUBLE -c -o $@ $^

%-CPU: $(CPUOBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ *-CPU.o

%-MIC: $(MICOBJECTS)
	$(CXX) $(CXXFLAGS) -mmic -o $@ *-MIC.o

all: $(TARGETS)

clean: 
	rm -f $(TARGETS) $(CPUOBJECTS) $(MICOBJECTS)
