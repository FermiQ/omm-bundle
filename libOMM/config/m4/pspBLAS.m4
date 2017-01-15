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
# pspBLAS support
#



# OMM_PSP_DETECT()
# -------------------
#
# Checks that the selected linear algebra libraries properly work.
#
AC_DEFUN([OMM_PSP_DETECT],[
  dnl Init
  omm_psp_ok="unknown"

  dnl Prepare environment
  saved_CPPFLAGS="${CPPFLAGS}"
  saved_FCFLAGS="${FCFLAGS}"
  saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${omm_psp_incs}"
  FCFLAGS="${FCFLAGS} ${CPPFLAGS}"
  LIBS="${omm_psp_libs} ${LIBS}"

  dnl Check pspBLAS routine
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether the pspBLAS library works])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use pspblas
      integer :: M,N,K
      type(psp_matrix_spm) :: A
      character(1) :: opA, opB
      double precision :: B(2,2),C(2,2),alpha,beta
      call psp_gespmm(M,N,K,A,opA,B,opB,C,alpha,beta)
    ]])], [omm_psp_ok="yes"], [omm_psp_ok="no"])
  AC_MSG_RESULT([${omm_psp_ok}])
  AC_LANG_POP([Fortran])

  dnl Restore environment
  CPPFLAGS="${saved_CPPFLAGS}"
  FCFLAGS="${saved_FCFLAGS}"
  LIBS="${saved_LIBS}"
]) # OMM_PSP_DETECT
