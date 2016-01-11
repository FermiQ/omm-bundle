PSPDIR = ./
include make.inc

.SUFFIXES:
.SUFFIXES: .o .a .f90 .F90

all : clean update pspBLASMAIN $(LIBNAME)

pspBLASMAIN : pspBLAS.o; \
              mv *.o ./include; \
              mv *.mod ./include

ALLOBJ = ./include/*.o

$(LIBNAME) : 
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@; \
        mv *.a ./include

update :
	cd ./pspMPI; \
	make update; \
	cd ../pspVariable; \
	make update; \
	cd ../pspUtility; \
	make update; \
	cd ../pspLevel1; \
	make update; \
	cd ../pspLevel2; \
	make update; \
	cd ../pspLevel3; \
	make update; \
	cd ..

clean : 
	rm -f *.a; \
	rm -f *.o; \
	rm -f *.mod; \
	cd ./include; \
	rm -f *.o; \
	rm -f *.mod; \
	rm -f *.a; \
	cd ../test/; \
	rm -f *.x

.F90.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $< -I$(PSPINCLUDEPATH) -o $@

.f90.o : 
	$(FORTRAN) $(OPTS) -c $< -I$(PSPINCLUDEPATH) -o $@

