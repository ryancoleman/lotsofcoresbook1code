include make.inc

FSRC = \
	   scalapack_mod.F \
	   crandom_number.F \
	   assert.F \
	   prof_mod.F \
	   pgesvr_mod.F 

FOBJ = \
	   scalapack_mod.o \
	   crandom_number.o \
	   assert.o \
	   prof_mod.o \
	   pgesvr_mod.o 



ISRC = \
	Cindxl2g.cpp \
	Cindxg2p.cpp \
	Cinfog1l.cpp \
	Cinfog2l.cpp \
	Cdescinit.cpp \
	Cdescset.cpp \
	Cnumroc2.cpp \
	Cnumroc.cpp \
	local_extent.cpp \
	Ciafirst.cpp \
	Cialast.cpp \
	Cpilaprnt.cpp \
	profinit.c \
	setup_desc.cpp \
	MallocHost.cpp \
	cublasinit.cpp \

IOBJ = \
	Cindxl2g.o \
	Cindxg2p.o \
	Cinfog1l.o \
	Cinfog2l.o \
	Cdescinit.o \
	Cdescset.o \
	Cnumroc2.o \
	Cnumroc.o \
	local_extent.o \
	Ciafirst.o \
	Cialast.o \
	Cpilaprnt.o \
	profinit.o \
	setup_desc.o \
	MallocHost.o \
	cublasinit.o \

ZSRC = \
	Cpzlaprnt.cpp \
	Cpzgecopy_d2h.cpp \
	Cpzgecopy_d2h_async.cpp \
	Cpzgemm_hhd.cpp \
	Cpzlaswp_gpu.cpp \
	Cpzswap_gpu.cpp \
	pzgetrf_gpu.cpp \
	pzgetrf_gpu2.cpp \
	pzgetrf_ooc.cpp \
	pzgetrf_ooc2.cpp 


ZOBJ = \
	Cpzlaprnt.o \
	Cpzgecopy_d2h.o \
	Cpzgecopy_d2h_async.o \
	Cpzgemm_hhd.o \
	Cpzlaswp_gpu.o \
	Cpzswap_gpu.o \
	pzgetrf_gpu.o \
	pzgetrf_gpu2.o \
	pzgetrf_ooc.o \
	pzgetrf_ooc2.o 

CSRC = \
	Cpclaprnt.cpp \
	Cpcgecopy_d2h.cpp \
	Cpcgecopy_d2h_async.cpp \
	Cpcgemm_hhd.cpp \
	Cpclaswp_gpu.cpp \
	Cpcswap_gpu.cpp \
	pcgetrf_gpu.cpp \
	pcgetrf_gpu2.cpp \
	pcgetrf_ooc.cpp \
	pcgetrf_ooc2.cpp 



COBJ = \
	Cpclaprnt.o \
	Cpcgecopy_d2h.o \
	Cpcgecopy_d2h_async.o \
	Cpcgemm_hhd.o \
	Cpclaswp_gpu.o \
	Cpcswap_gpu.o \
	pcgetrf_gpu.o \
	pcgetrf_gpu2.o \
	pcgetrf_ooc.o \
	pcgetrf_ooc2.o 


DSRC = \
	Cpdlaprnt.cpp \
	Cpdgecopy_d2h.cpp \
	Cpdgecopy_d2h_async.cpp \
	Cpdgemm_hhd.cpp \
	Cpdsyrk_hhd.cpp \
	pdpotrf_gpu2.cpp \
	pdpotrf_ooc2.cpp \
#	Cpdswap_gpu.cpp \
#	Cpdlaswp_gpu.cpp \
#	pdcopymatrix.cpp \
#	pdgetrf_gpu.cpp \
#	pdgetrf_gpu2.cpp \
#	pdgetrf_ooc.cpp \
#	pdgetrf_ooc2.cpp 
        


DOBJ = \
	Cpdlaprnt.o \
	Cpdgecopy_d2h.o \
	Cpdgecopy_d2h_async.o \
	Cpdgemm_hhd.o \
	Cpdsyrk_hhd.o \
	pdpotrf_gpu2.o \
	pdpotrf_ooc2.o \
#	Cpdlaswp_gpu.o \
#	Cpdswap_gpu.o \
#	pdcopymatrix.o \
#	pdgetrf_gpu.o \
#	pdgetrf_gpu2.o \
#	pdgetrf_ooc.o \
#	pdgetrf_ooc2.o  

SSRC = \
	Cpslaprnt.cpp \
	Cpsgecopy_d2h.cpp \
	Cpsgecopy_d2h_async.cpp \
	Cpsgemm_hhd.cpp \
	Cpssyrk_hhd.cpp \
	Cpslaswp_gpu.cpp \
	Cpsswap_gpu.cpp \
	psgetrf_gpu.cpp \
	psgetrf_gpu2.cpp \
	psgetrf_ooc.cpp \
	psgetrf_ooc2.cpp \
        pspotrf_ooc2.cpp \
        pspotrf_gpu2.cpp



SOBJ = \
	Cpslaprnt.o \
	Cpsgecopy_d2h.o \
	Cpsgecopy_d2h_async.o \
	Cpsgemm_hhd.o \
	Cpssyrk_hhd.o \
	Cpslaswp_gpu.o \
	Cpsswap_gpu.o \
	psgetrf_gpu.o \
	psgetrf_gpu2.o \
	psgetrf_ooc.o \
	psgetrf_ooc2.o \
	pspotrf_ooc2.o \
	pspotrf_gpu2.o


SRC= $(ISRC) $(DSRC)#$(CSRC) $(FSRC) $(ZSRC) $(SSRC)

OBJECTS = $(IOBJ) $(DOBJ)#$(COBJ) $(FOBJ) $(ZOBJ) $(SOBJ)


include Make.lindep
liblin.a: $(LINSRC)
	$(FC) $(FFLAGS) -c $(LINSRC)
	ar r liblin.a $(LINOBJECTS); ranlib liblin.a

lib.a: $(OBJECTS)

	 ar r lib.a $(OBJECTS); ranlib lib.a
	  
profinit.o: profinit.c
	$(CC) $(CFLAGS) -c profinit.c

Cindxl2g.o: Cindxl2g.cpp
	$(CC) $(CFLAGS) -c Cindxl2g.cpp

Cindxg2p.o: Cindxg2p.cpp
	$(CC) $(CFLAGS) -c Cindxg2p.cpp

Cdescset.o: Cdescset.cpp
	$(CC) $(CFLAGS) -c Cdescset.cpp

Cinfog2l.o: Cinfog2l.cpp
	$(CC) $(CFLAGS) -c Cinfog2l.cpp


Cialast.o: Cialast.cpp
	$(CC) $(CFLAGS) -c Cialast.cpp

Ciafirst.o: Ciafirst.cpp
	$(CC) $(CFLAGS) -c Ciafirst.cpp

local_extent.o: local_extent.cpp
	$(CC) $(CFLAGS) -c local_extent.cpp

MallocHost.o: MallocHost.cpp
	$(CC) $(CFLAGS) -c MallocHost.cpp

setup_desc.o: setup_desc.cpp
	$(CC) $(CFLAGS) -c setup_desc.cpp

Cdescinit.o: Cdescinit.cpp
	$(CC) $(CFLAGS) -c Cdescinit.cpp


Cinfog1l.o: Cinfog1l.cpp
	$(CC) $(CFLAGS) -c Cinfog1l.cpp

Cpilaprnt.o: Cpilaprnt.cpp
	$(CC) $(CFLAGS) -c Cpilaprnt.cpp


Cnumroc.o: Cnumroc.cpp
	$(CC) $(CFLAGS) -c Cnumroc.cpp

Cnumroc2.o: Cnumroc2.cpp
	$(CC) $(CFLAGS) -c Cnumroc2.cpp

cublasinit.o: cublasinit.cpp
	$(CC) $(CFLAGS) -c cublasinit.cpp

test_Cialast: lib.a test_Cialast.cpp Cialast.cpp
	$(CC) $(CFLAGS) -o test_Cialast \
	  	test_Ciafirst.cpp \
		lib.a $(LIBS)

test_Ciafirst: lib.a test_Ciafirst.cpp Ciafirst.cpp
	$(CC) $(CFLAGS) -o test_Ciafirst \
	  	test_Ciafirst.cpp \
		lib.a $(LIBS)

crandom_number.o: crandom_number.F
	$(FC) $(FFLAGS) -c crandom_number.F

assert.o: assert.F
	$(FC) $(FFLAGS) -c assert.F

prof_mod.o: prof_mod.F
	$(FC) $(FFLAGS) -c prof_mod.F


scalapack_mod.o: scalapack_mod.F
	$(FC) $(FFLAGS) -c scalapack_mod.F

pgesvr_mod.o: scalapack_mod.F pgesvr_mod.F 
	$(FC) $(FFLAGS) -c pgesvr_mod.F


tpdgesvr.exe: $(FOBJ) lib.a tpdgesvr.F
	$(FC) $(FFLAGS) -o tpdgesvr.exe \
	  tpdgesvr.F \
	  lib.a $(LIBS)

tpzgesvr.exe: $(FOBJ) lib.a tpzgesvr.F
	$(FC) $(FFLAGS) -o tpzgesvr.exe \
	  tpzgesvr.F \
	  lib.a $(LIBS)

include makedepend_z
include makedepend_c
include makedepend_d
include makedepend_s
clean:
	touch fake_cublas.o  lib.a
	touch pdludriver.exe pdludriver_fake.exe  pdludriver2.exe
	touch pzludriver.exe pzludriver_fake.exe  pzludriver2.exe
	touch pcludriver.exe pcludriver_fake.exe  pcludriver2.exe
	touch test_Ciafirst test_Cialast
	rm -f *.o  lib.a 
	rm -f pdllt*.exe
	rm -f pdlu*.exe
	rm -f pzlu*.exe
	rm -f pclu*.exe
	rm -f test_Ciafirst test_Cialast
	rm -f liblin.a
