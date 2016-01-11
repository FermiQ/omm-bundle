include make.inc

.SUFFIXES:
.SUFFIXES: .x .f90 .F90

all : pdgemmScaling.x pzgemmScaling.x

pdgemmScaling.x : pdgemmScaling.F90

pzgemmScaling.x : pzgemmScaling.F90

update :
	cd ../src/; \
	make clean; \
	make; \
	cd ../examples; \
	make clean; \
	make; \
	mpirun -n 4 ./pzgemmScaling.x

clean : 
	rm -f *.x

LINK_MACRO = $< -I../include ../include/libpspblas.a $(SLAPLIB) $(LAPLIB) $(BLASLIB) -o $@

.F90.x : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) $(LINK_MACRO)

.f90.x : 
	$(FORTRAN) $(OPTS) $(LINK_MACRO)
