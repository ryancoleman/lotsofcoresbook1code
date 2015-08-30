.SUFFIXES: .c .o

DIR = ~/nbody/mic
LIB = /opt/intel/composer_xe_2013.1.117/compiler/lib/mic
MIC = mic0

all:
	make clean
	make build
	make run

build: nbody.c
	icc -mmic -openmp -fimf-domain-exclusion=15 -O3 $<
run:
	scp $(LIB)/libiomp5.so $(MIC):$(DIR)
	ssh $(MIC) LD_LIBRARY_PATH=$(DIR) $(DIR)/a.out
clean:
	find . -name "*.o" -o -name "*.out*" | xargs rm -rf
