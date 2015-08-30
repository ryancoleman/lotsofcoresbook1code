CXXFLAGS=-mkl -openmp -fimf-precision=low -opt-report -opt-report-phase:vec -opt-report-stdout -O3
MICFLAGS=-mmic -fimf-domain-exclusion=15
HOSTFLAGS=-xHOST

all:
	mpiicpc $(CXXFLAGS) $(MICFLAGS)  -o options.MIC options.c
	mpiicpc $(CXXFLAGS) $(HOSTFLAGS) -o options     options.c

run: run-cpu run-mic run-both

run-cpu:
	./run-options.sh cpu

run-mic:
	./run-options.sh mic

run-both:
	./run-options.sh both

clean:
	rm -f options options.MIC results* *.stf options.optrpt  options.prot
