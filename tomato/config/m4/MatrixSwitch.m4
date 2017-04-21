# -*- Autoconf -*-
#
# M4 macros for Tomato
#
# Copyright (C) 2017 Yann Pouillon
#
# This file is part of the Tomato software package. For license information,
# please see the COPYING file in the top-level directory of the Tomato source
# distribution.
#

#
# MatrixSwitch support
#



# TOM_MSW_DETECT()
# -------------------
#
# Checks that the selected linear algebra libraries properly work.
#
AC_DEFUN([TOM_MSW_DETECT],[
  dnl Init
  tom_msw_ok="unknown"

  dnl Prepare environment
  saved_CPPFLAGS="${CPPFLAGS}"
  saved_FCFLAGS="${FCFLAGS}"
  saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${tom_msw_incs}"
  FCFLAGS="${FCFLAGS} ${CPPFLAGS}"
  LIBS="${tom_msw_libs} ${LIBS}"

  dnl Check MatrixSwitch routine
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether the MatrixSwitch library works])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use matrixswitch
      character(3) :: m_operation
      double precision :: x_min
      type(matrix) :: HW, HWd
      call m_add(HWd, 'c', HW, x_min, 1.0d0, m_operation)
    ]])], [tom_msw_ok="yes"], [tom_msw_ok="no"])
  AC_MSG_RESULT([${tom_msw_ok}])
  AC_LANG_POP([Fortran])

  dnl Restore environment
  CPPFLAGS="${saved_CPPFLAGS}"
  FCFLAGS="${saved_FCFLAGS}"
  LIBS="${saved_LIBS}"
]) # TOM_MSW_DETECT
