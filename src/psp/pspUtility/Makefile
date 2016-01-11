PSPDIR = ../
include $(PSPDIR)/make.inc

.SUFFIXES:
.SUFFIXES: .o .a .f90 .F90

all : pspUtility.o 

clean : 
	rm -f *.a
	rm -f *.o
	rm -f *.mod

update :
	make clean; \
	make; \
	mv *.o ../include; \
	mv *.mod ../include


.F90.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $< -I$(PSPINCLUDEPATH) -o $@

.f90.o : 
	$(FORTRAN) $(OPTS) -c $< -I$(PSPINCLUDEPATH) -o $@
