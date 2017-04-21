# -*- Autoconf -*-
#
# M4 macros for Tomato
#
# Copyright (C) 2016 Yann Pouillon
#
# This file is part of the Tomato software package. For license information,
# please see the COPYING file in the top-level directory of the Tomato source
# distribution.
#

#
# Linear algebra support
#



# TOM_LINALG_DETECT()
# -------------------
#
# Checks that the selected linear algebra libraries properly work.
#
AC_DEFUN([TOM_LINALG_DETECT],[
  dnl Init
  tom_linalg_has_lapack="unknown"
  tom_linalg_has_scalapack="unknown"
  tom_linalg_ok="unknown"

  dnl Prepare environment
  saved_CPPFLAGS="${CPPFLAGS}"
  saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${tom_linalg_incs}"
  LIBS="${tom_linalg_libs} ${LIBS}"

  dnl Check BLAS and LAPACK routines
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether linear algebra libraries work])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zgemm
      call zhpev
    ]])], [tom_linalg_ok="yes"; tom_linalg_has_lapack="yes"], [tom_linalg_ok="no"])
  AC_MSG_RESULT([${tom_linalg_ok}])
  AC_LANG_POP([Fortran])

  dnl Check ScaLAPACK routines
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether linear algebra libraries have ScaLAPACK])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call pzheevx
    ]])], [tom_linalg_has_scalapack="yes"], [tom_linalg_has_scalapack="no"])
  AC_MSG_RESULT([${tom_linalg_ok}])
  AC_LANG_POP([Fortran])

  dnl Restore environment
  CPPFLAGS="${saved_CPPFLAGS}"
  LIBS="${saved_LIBS}"
]) # TOM_LINALG_DETECT
