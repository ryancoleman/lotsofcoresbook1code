#
#  USAGE:
#     make          ... to build the programs
#     make clean    ... to remove object and executable files
#

# verify that you are including the right make.def file for the platform
include make.def

EXES     = matmul$(EXE)

MMULOBJS = matmul_host.$(OBJ) matmul_lib.$(OBJ)

all: $(EXES)

matmul$(EXE): $(MMULOBJS) 
	$(CLINKER) $(CFLAGS) $(OPENCLFLAGS) -o $@ $^ $(LIBS)

$(MMULOBJS): matmul.hpp cl.hpp util.hpp

matmul_host.$(OBJ): matmul_lib.hpp

clean:
	$(RM) $(EXES) *.$(OBJ)

veryclean:
	$(RM) $(EXES) *.$(OBJ)

.SUFFIXES:
.SUFFIXES: .c .cpp .$(OBJ)

.c.$(OBJ):
	$(CC) $(CFLAGS) -c $<

.cpp.$(OBJ):
	$(CC) $(CFLAGS) -c $<


