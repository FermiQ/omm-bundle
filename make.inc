SHELL     = /bin/sh
FORTRAN   = mpif90
OPTS      = -O3
FPPFLAGS  = -DCONV -DPSP -DMPI -DLAP -DSLAP -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

PSPLIBNAME = libpspblas.a
PSPBUILDPATH = $(ROOTDIR)/build_psp

MSLIBNAME = MatrixSwitch.a
MSBUILDPATH = $(ROOTDIR)/build_MatrixSwitch

OMMLIBNAME = libOMM.a
OMMBUILDPATH = $(ROOTDIR)/build_libOMM

TOMATOLIBNAME = tomato.a
TOMATOBUILDPATH = $(ROOTDIR)/build_tomato

BLASLIB = 
LAPLIB  = 
SLAPLIB = ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl
