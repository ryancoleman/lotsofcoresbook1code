CC=icpc
FLAGS=-O3 -fp-model fast=2 -mmic -openmp -mkl=sequential
MIC=mic0

BIN=nbody

all:$(BIN) run

nbody: nbody.cc
	$(CC) -o $@  $(FLAGS) $<

run:
	@ssh $(MIC) LD_LIBRARY_PATH=$$MIC_LD_LIBRARY_PATH $$(pwd)/nbody

clean:
	rm -f $(BIN)
