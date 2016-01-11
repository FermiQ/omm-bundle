SHELL     = /bin/sh
FORTRAN   = mpif90
OPTS      = -O3
ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib
FPPFLAGS := $(FPPFLAGS) -DCONV -DMPI -DLAP -DSLAP
LIBNAME   = libpspblas.a

PSPINCLUDEPATH = $(PSPDIR)/include/
