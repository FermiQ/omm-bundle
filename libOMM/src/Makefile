include make.inc

.SUFFIXES:
.SUFFIXES: .o .a .f90 .F90

all : $(LIBNAME)

ALLOBJ = omm_params.o omm_rand.o omm_quartic.o omm_ops.o omm.o omm_wrappers.o omm_wrappers_spm.o

$(LIBNAME) : $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

install :
	if [ ! -d $(BUILDPATH) ]; then mkdir $(BUILDPATH); fi; \
	if [ ! -d $(BUILDPATH)/include ]; then mkdir $(BUILDPATH)/include; fi; \
	if [ ! -d $(BUILDPATH)/lib ]; then mkdir $(BUILDPATH)/lib; fi; \
	cp *.mod $(BUILDPATH)/include; \
	cp *.a $(BUILDPATH)/lib;

clean : 
	rm -f *.a
	rm -f *.o
	rm -f *.mod

.F90.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $< -I$(MSINCLUDEPATH) -I$(PSPINCLUDEPATH) -o $@

.f90.o : 
	$(FORTRAN) $(OPTS) -c $< -I$(MSINCLUDEPATH) -I$(PSPINCLUDEPATH) -o $@
