LIBPATH   = $(BASEDIR)
MSLIBPATH = $(BASEDIR)/MatrixSwitch-0.1.2/src/
PSPINCLUDE = $(BASEDIR)/psp/include

MSLIB     = $(BASEDIR)/MatrixSwitch-0.1.2/src/MatrixSwitch.a
PSPLIB   = $(PSPINCLUDE)/libpspblas.a
BLASLIB   = -lblas
LAPLIB    = -llapack
SLAPLIB   = $(BASEDIR)/../../scalapack/libscalapack.a