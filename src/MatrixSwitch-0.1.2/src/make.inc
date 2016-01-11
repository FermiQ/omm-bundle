include ../../make.inc

FPPFLAGS := $(FPPFLAGS) -DCONV
LIBNAME   = MatrixSwitch.a
PSPINCLUDE = ../../psp/include
PSPLIB   = $(PSPINCLUDE)/libpspblas.a
