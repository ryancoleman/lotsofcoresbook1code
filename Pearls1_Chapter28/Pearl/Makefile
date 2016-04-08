
all: miczmatmul zmatmul miczmattr zmattr

zmattr: zmattr.c zorder2d.c timer.c
	icc -O3 -fp-model fast=2 -openmp -lrt $^ -o $@

miczmattr: zmattr.c zorder2d.c timer.c
	icc -mmic -O3 -fp-model fast=2 -openmp -lrt $^ -o $@

zmatmul: zmatmul.c zorder2d.c timer.c
	icc -O3 -fp-model fast=2 -openmp -lrt $^ -o $@

miczmatmul: zmatmul.c zorder2d.c timer.c
	icc -mmic -O3 -fp-model fast=2 -openmp -lrt $^ -o $@

clean:
	rm zmatmul zmattr miczmatmul miczmattr
