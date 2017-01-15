omm-bundle
==========

This repository is a bundle of four separate libraries:

*   pspBLAS: MPI parallel sparse matrix operations with 2D block-cyclic
             distribution
*   MatrixSwitch: multi-format matrix storage and operation
*   libOMM: Kohn-Sham solver using the orbital minimization method
*   tomato: test matrix generator for Kohn-Sham solvers

Each library has its own documentation and build systems, the details of which
you can find in the individual subdirectories.

Additionally, in this directory you will find two ways of installing the library
bundle as a whole. They are explained below.

Autotools build system
----------------------

The Autotools build system is available for pspBLAS, MatrixSwitch, and libOMM
(tomato coming soon). It can be accessed by running the `build-omm` script and
following the instructions. Note that external linear algebra libraries need to
be given either by modifying the script or by creating environmental variables
by typing:

`export LINALG_INCLUDES="path/to/headers"`
`export LINALG_LIBS="path/to/libraries"`

Manual build system
-------------------

The manual build system is available for all libraries, and makes use of the
file named Makefile.manual.

To install the libraries manually:

1.  Copy `make.inc.example` to `make.inc` and modify it to suit your needs.
    Available options for `FPPFLAGS` are:
    * `-DHAVE_MPI`: enable MPI parallel routines
    * `-DHAVE_LAPACK`: enable LAPACK routines
    * `-DHAVE_SCALAPACK`: enable ScaLAPACK routines (requires MPI)
    * `-DHAVE_PSPBLAS`: enable pspBLAS routines (requires MPI, LAPACK and
      ScaLAPACK)
    * `-DCONV`: enable automatic conversion of scalar types (real/complex) to
      agree with matrix definitions (real/complex). Needed for libOMM.
    * `-DNORAND`: fixed seed for the random number generator. Enable for testing
      purposes.
    * `-DCBIND`: use ISO_C_BINDING for LOGICAL inputs in the wrapper interfaces.
      Enable for linking to C.
3.  Type `make -f Makefile.manual`.
