.SUFFIXES: .c .o

all:
	make clean
	make build
	make run

build: nbody.c
	icc -mavx -openmp -O3 $<
run:
	./a.out
clean:
	find . -name "*.o" -o -name "*.out*" | xargs rm -rf