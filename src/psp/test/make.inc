SHELL     = /bin/sh
FORTRAN   = mpif90
OPTS      = -O3
FPPFLAGS  = -DMPI -DLAP -DSLAP

BLASLIB   = -lblas
LAPLIB    = -llapack
SLAPLIB   = ../../../../scalapack/libscalapack.a