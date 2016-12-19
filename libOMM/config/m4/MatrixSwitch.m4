# -*- Autoconf -*-
#
# M4 macros for libOMM
#
# Copyright (C) 2016 Yann Pouillon
#
# This file is part of the libOMM software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

#
# MatrixSwitch support
#



# OMM_MSW_DETECT()
# -------------------
#
# Checks that the selected linear algebra libraries properly work.
#
AC_DEFUN([OMM_MSW_DETECT],[
  dnl Init
  omm_matrixswitch_ok="unknown"

  dnl Prepare environment
  saved_CPPFLAGS="${CPPFLAGS}"
  saved_FCFLAGS="${FCFLAGS}"
  saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${omm_matrixswitch_incs}"
  FCFLAGS="${FCFLAGS} ${omm_matrixswitch_incs}"
  LIBS="${omm_matrixswitch_libs} ${LIBS}"

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
    ]])], [omm_matrixswitch_ok="yes"], [omm_matrixswitch_ok="no"])
  AC_MSG_RESULT([${omm_matrixswitch_ok}])
  AC_LANG_POP([Fortran])

  dnl Restore environment
  CPPFLAGS="${saved_CPPFLAGS}"
  FCFLAGS="${saved_FCFLAGS}"
  LIBS="${saved_LIBS}"
]) # OMM_MSW_DETECT
