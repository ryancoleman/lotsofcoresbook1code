#
# $Id: GNUmakefile 19707 2010-10-29 17:59:36Z d3y133 $
#


           OBJ = md5wrap.o md5.o
       LIBRARY = libnwcutil.a
   LIB_TARGETS = ftest mddriver
       HEADERS = checksum.h

include ../../config/makefile.h
include ../../config/makelib.h

            COPTIMIZE += -O0 -g
            COPTIONS += -O0 -g

ftest:	ftest.o $(OBJ)
	$(LINK.f) -o $@ $^ $(LIBS) $(LIBS) $(LIBS) -lm

mddriver:	mddriver.o $(OBJ)
	$(CC) -o $@ $^ -lm

