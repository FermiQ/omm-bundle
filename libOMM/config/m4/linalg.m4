# -*- Autoconf -*-
#
# M4 macros for libOMM
#
# Copyright (C) 2016 Yann Pouillon
#
# This file is part of the libOMM software package. For license information,
# please see the COPYING file in the top-level directory of the libOMM source
# distribution.
#

#
# Linear algebra support
#



# OMM_LINALG_DETECT()
# -------------------
#
# Checks that the selected linear algebra libraries properly work.
#
AC_DEFUN([OMM_LINALG_DETECT],[
  dnl Init
  omm_linalg_has_lapack="unknown"
  omm_linalg_has_scalapack="unknown"
  omm_linalg_ok="unknown"

  dnl Prepare environment
  saved_CPPFLAGS="${CPPFLAGS}"
  saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${omm_linalg_incs}"
  LIBS="${omm_linalg_libs} ${LIBS}"

  dnl Check BLAS and LAPACK routines
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether linear algebra libraries work])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zgemm
      call zhpev
    ]])], [omm_linalg_ok="yes"; omm_linalg_has_lapack="yes"], [omm_linalg_ok="no"])
  AC_MSG_RESULT([${omm_linalg_ok}])
  AC_LANG_POP([Fortran])

  dnl Check ScaLAPACK routines
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether linear algebra libraries have ScaLAPACK])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call pzheevx
    ]])], [omm_linalg_has_scalapack="yes"], [omm_linalg_has_scalapack="no"])
  AC_MSG_RESULT([${omm_linalg_ok}])
  AC_LANG_POP([Fortran])

  dnl Restore environment
  CPPFLAGS="${saved_CPPFLAGS}"
  LIBS="${saved_LIBS}"
]) # OMM_LINALG_DETECT
