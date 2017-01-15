# Copyright (C) 2016 Yann Pouillon <notifications@materialsevolution.es>
#
# This file is part of the libOMM software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

#
# Configure script for the libOMM package
#

                    # ------------------------------------ #

#
# Autotools startup
#

# Init Autoconf
AC_PREREQ(2.69)
AC_INIT([libOMM], [0.0.0], [http://www.e-cam2020.eu:10080/ESL/omm/issues])
AC_REVISION([Autotools support for libOMM])
AC_CONFIG_AUX_DIR(config/gnu)
AC_CONFIG_MACRO_DIR([config/m4])
AC_CONFIG_SRCDIR([src/omm.F90])
_AC_SRCDIRS([.])

# Init Automake
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([1.14 parallel-tests color-tests])
AM_CONFIG_HEADER([config.h])

# Generate timestamp
omm_timestamp=`date '+%Y%m%dT%H%M%S%z'`
AC_SUBST(omm_timestamp)

                    # ------------------------------------ #

#
# System utilities
#

# Check for common programs
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_SED
AC_PROG_AWK
AC_PROG_GREP

                    # ------------------------------------ #

#
# Default settings
#

# Note for developers: you may edit these parameters to tune the behaviour
# of the build system. Please do it with care.

# Lists of options to manage (omm_*_enable_def variables)
#
# Note:
#   * switches (--enable-*) must be set to "yes" or "no"
#   * external dependencies (--with-*) must be set to "yes", "no", or "auto"
#
omm_cfg_switches="debug examples randseed"
omm_cfg_extdeps="linalg mpi psp msw"

# Optional features
omm_debug_enable_def="yes"
omm_examples_enable_def="no"
omm_randseed_enable_def="yes"

# Linear algebra
# Note: default libraries must be empty because some compilers have
#       built-in support.
omm_linalg_enable_def="no"
omm_linalg_incs_def=""
omm_linalg_libs_def=""

# MPI
omm_mpi_enable_def="auto"
omm_mpi_incs_def=""
omm_mpi_libs_def=""

# pspBLAS
omm_psp_enable_def="no"
omm_psp_incs_def=""
omm_psp_libs_def=""

# MatrixSwitch
omm_msw_enable_def="yes"
omm_msw_incs_def=""
omm_msw_libs_def=""

                    # ------------------------------------ #

#
# Check consistency of defaults
#

# Default triggers for optional features must be yes or no, and not empty
for omm_cfg_default in ${omm_cfg_switches}; do
  tmp_default=`eval echo \$\{omm_${omm_cfg_default}_enable_def\}`
  AC_MSG_CHECKING([default value of the ${omm_cfg_default} switch])
  AC_MSG_RESULT([${tmp_default}])
  if test "${tmp_default}" != "no" -a \
          "${tmp_default}" != "yes"; then
    AC_MSG_ERROR([invalid default value: omm_${omm_cfg_default}_enable_def = ${tmp_default}])
  fi
done

# Default triggers for external dependencies must be yes, no, or auto
for omm_cfg_default in ${omm_cfg_extdeps}; do
  tmp_default=`eval echo \$\{omm_${omm_cfg_default}_enable_def\}`
  AC_MSG_CHECKING([default detection behavior for ${omm_cfg_default}])
  AC_MSG_RESULT([${tmp_default}])
  if test "${tmp_default}" != "" -a \
          "${tmp_default}" != "auto" -a \
          "${tmp_default}" != "no" -a \
          "${tmp_default}" != "yes"; then
    AC_MSG_ERROR([invalid default value: omm_${omm_cfg_default}_enable_def = ${tmp_default}])
  fi
done

# Clean-up
unset tmp_default

                    # ------------------------------------ #

#
# Optional features
#

# Conditional compilation: --enable-* options (alphabetical order)
AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug],
    [Enable extra debug information])],
  [omm_debug_enable="${enableval}"; omm_debug_type="yon"],
  [omm_debug_enable="${omm_debug_enable_def}"; omm_debug_type="def"])
AC_SUBST(enable_debug)
AC_ARG_ENABLE([examples],
  [AS_HELP_STRING([--enable-examples],
    [Build usage examples of libOMM])],
  [omm_examples_enable="${enableval}"; omm_examples_type="yon"],
  [omm_examples_enable="${omm_examples_enable_def}"; omm_examples_type="def"])
AC_SUBST(enable_examples)
AC_ARG_ENABLE([randseed],
  [AS_HELP_STRING([--enable-randseed],
    [Enable time-based random seed, or else 123456])],
  [omm_randseed_enable="${enableval}"; omm_randseed_type="yon"],
  [omm_randseed_enable="${omm_randseed_enable_def}"; omm_randseed_type="def"])
AC_SUBST(enable_randseed)

                    # ------------------------------------ #

#
# External dependencies
#

# MPI (optional)
AC_ARG_WITH([mpi],
  [AS_HELP_STRING([--with-mpi],
    [Install prefix of MPI (e.g. /usr/local). The default behaviour is to detect whether the specified compilers properly support MPI and to fall back to serial mode if not. You may use --with-mpi without argument to force MPI detection, in which case detection failures will result in errors, and --without-mpi to disable MPI support completely.])],
  [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
      omm_mpi_enable="${withval}"
      omm_mpi_type="yon"
    else
      omm_mpi_enable="yes"
      omm_mpi_type="dir"
    fi],
  [omm_mpi_enable="${omm_mpi_enable_def}"; omm_mpi_type="def"])
AC_SUBST(with_mpi)

AC_ARG_VAR([MPICC], [MPI-aware C compiler (usually called 'mpicc'). Please note that this variable conflicts with the --with-mpi option.])
if test -n "${MPICC}"; then
  omm_mpi_enable="yes"
  omm_mpi_type="env"
fi
AC_ARG_VAR([MPIFC], [MPI-aware Fortran compiler (usually called 'mpif90'. Please note that this variable conflicts with the --with-mpi option.])
if test -n "${MPIFC}"; then
  omm_mpi_enable="yes"
  omm_mpi_type="env"
fi
AC_ARG_VAR([MPI_INCLUDES],
  [Include flags for MPI. Please note that this variable conflicts with the --with-mpi option.])
if test -n "${MPI_INCLUDES}"; then
  omm_mpi_enable="yes"
  omm_mpi_type="env"
  omm_mpi_incs="${MPI_INCLUDES}"
fi
AC_ARG_VAR([MPI_LIBS],
    [Library flags for MPI. Please note that this variable conflicts with the --with-mpi option.])
if test -n "${MPI_LIBS}"; then
  omm_mpi_enable="yes"
  omm_mpi_type="env"
  omm_mpi_libs="${MPI_LIBS}"
fi

# Linear algebra (optional)
AC_ARG_WITH([linalg],
  [AS_HELP_STRING([--with-linalg],
    [Install prefix of the linear algebra libraries (e.g. /usr/local). Linear algebra support will be enabled only if you use this option. You may use --with-linalg without argument if linear algebra libraries are available system-wide, and --without-linalg to disable linear algebra support completely.])],
  [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
      omm_linalg_enable="${withval}"
      omm_linalg_type="yon"
    else
      omm_linalg_enable="yes"
      omm_linalg_type="dir"
      omm_linalg_incs="-I${withval}/include"
      omm_linalg_libs="-L${withval}/lib"
    fi],
  [omm_linalg_enable="${omm_linalg_enable_def}"; omm_linalg_type="def";
    omm_linalg_incs=""; omm_linalg_libs=""])
AC_SUBST(with_linalg)

AC_ARG_VAR([LINALG_INCLUDES],
  [Include flags for linear algebra. Please note that this variable conflicts with the --with-linalg option.])
if test -n "${LINALG_INCLUDES}"; then
  omm_linalg_enable="yes"
  omm_linalg_type="env"
fi
AC_ARG_VAR([LINALG_LIBS],
    [Library flags for linear algebra. Please note that this variable conflicts with the --with-linalg option.])
if test -n "${LINALG_LIBS}"; then
  omm_linalg_enable="yes"
  omm_linalg_type="env"
fi

# pspBLAS (optional)
AC_ARG_WITH([psp],
  [AS_HELP_STRING([--with-psp],
    [Install prefix of the pspBLAS library (e.g. /usr/local). pspBLAS support will be enabled only if you use this option. You may use --with-psp without argument if pspBLAS is installed as a system library, and --without-psp to disable pspBLAS support completely.])],
  [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
      omm_psp_enable="${withval}"
      omm_psp_type="yon"
    else
      omm_psp_enable="yes"
      omm_psp_type="dir"
      omm_psp_incs="-I${withval}/include"
      omm_psp_libs="-L${withval}/lib -lpspBLAS"
    fi],
  [omm_psp_enable="${omm_psp_enable_def}"; omm_psp_type="def";
    omm_psp_incs=""; omm_psp_libs=""])
AC_SUBST(with_psp)

AC_ARG_VAR([PSP_INCLUDES],
  [Include flags for pspBLAS. Please note that this variable conflicts with the --with-psp option.])
if test -n "${PSP_INCLUDES}"; then
  omm_psp_enable="yes"
  omm_psp_type="env"
fi
AC_ARG_VAR([PSP_LIBS],
    [Library flags for pspBLAS. Please note that this variable conflicts with the --with-psp option.])
if test -n "${PSP_LIBS}"; then
  omm_psp_enable="yes"
  omm_psp_type="env"
fi

# MatrixSwitch (optional)
AC_ARG_WITH([msw],
  [AS_HELP_STRING([--with-msw],
    [Install prefix of the MatrixSwitch library (e.g. /usr/local). MatrixSwitch support will be enabled only if you use this option. You may use --with-msw without argument if MatrixSwitch is installed as a system library, and --without-msw to disable MatrixSwitch support completely.])],
  [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
      omm_msw_enable="${withval}"
      omm_msw_type="yon"
    else
      omm_msw_enable="yes"
      omm_msw_type="dir"
      omm_msw_incs="-I${withval}/include"
      omm_msw_libs="-L${withval}/lib -lMatrixSwitch"
    fi],
  [omm_msw_enable="${omm_msw_enable_def}"; omm_msw_type="def";
    omm_msw_incs=""; omm_msw_libs=""])
AC_SUBST(with_msw)

AC_ARG_VAR([MSW_INCLUDES],
  [Include flags for MatrixSwitch. Please note that this variable conflicts with the --with-msw option.])
if test -n "${MSW_INCLUDES}"; then
  omm_msw_enable="yes"
  omm_msw_type="env"
fi
AC_ARG_VAR([MSW_LIBS],
    [Library flags for MatrixSwitch. Please note that this variable conflicts with the --with-msw option.])
if test -n "${MSW_LIBS}"; then
  omm_msw_enable="yes"
  omm_msw_type="env"
fi

                    # ------------------------------------ #

#
# Check option consistency
#

# All --enable-* options must be yes or no
for omm_cfg_option in ${omm_cfg_switches}; do
  tmp_option=`eval echo \$\{enable_${omm_cfg_option}\}`
  if test "${tmp_option}" != "" -a \
          "${tmp_option}" != "no" -a \
          "${tmp_option}" != "yes"; then
    AC_MSG_ERROR([--enable-${omm_cfg_option} must be "yes" or "no"])
  fi
done

# Linear algebra must be enabled
#if test "${omm_linalg_enable}" != "yes"; then
#  AC_MSG_ERROR([linear algebra support has been disabled
#                  while it is mandatory])
#fi

# All MPI wrappers must be provided at once
tmp_mpi_vars="${MPICC}${MPIFC}"
if test -n "${tmp_mpi_vars}" -a \
        \( -z "${MPICC}" -o -z "${MPIFC}" \); then
  AC_MSG_ERROR([MPICC and MPIFC must both be set])
fi

# All pspBLAS environment variables must be provided at once
tmp_psp_vars="${PSP_INCLUDES}${PSP_LIBS}"
if test -n "${tmp_psp_vars}" -a \
        \( -z "${PSP_INCLUDES}" -o -z "${PSP_LIBS}" \); then
  AC_MSG_ERROR([PSP_INCLUDES and PSP_LIBS must both be set])
fi

# All MatrixSwitch environment variables must be provided at once
tmp_msw_vars="${MSW_INCLUDES}${MSW_LIBS}"
if test -n "${tmp_msw_vars}" -a \
        \( -z "${MSW_INCLUDES}" -o -z "${MSW_LIBS}" \); then
  AC_MSG_ERROR([MSW_INCLUDES and MSW_LIBS must both be set])
fi

# Environment variables conflict with --with-* options
# Note: for linear algebra, only LINALG_LIBS is relevant
if test -n "${tmp_mpi_vars}" -a -n "${with_mpi}"; then
  AC_MSG_ERROR([conflicting option settings
                  Please use MPICC + MPIFC or --with-mpi,
                  not both.])
fi
if test -n "${LINALG_LIBS}" -a -n "${with_linalg}"; then
  AC_MSG_ERROR([conflicting option settings
                  Please use LINALG_LIBS or --with-linalg,
                  not both.])
fi
if test -n "${tmp_psp_vars}" -a -n "${with_psp}"; then
  AC_MSG_ERROR([conflicting option settings
                  Please use PSP_INCLUDES + PSP_LIBS or --with-psp,
                  not both.])
fi

# When using environment variables, triggers must be set to yes
if test -n "${tmp_mpi_vars}"; then
  omm_mpi_enable="yes"
  omm_mpi_type="env"
fi
if test -n "${LINALG_LIBS}"; then
  omm_linalg_enable="yes"
  omm_linalg_type="env"
fi
if test -n "${tmp_psp_vars}"; then
  omm_psp_enable="yes"
  omm_psp_type="env"
fi
if test -n "${tmp_msw_vars}"; then
  omm_msw_enable="yes"
  omm_msw_type="env"
fi

# pspBLAS requires MPI
if test "${omm_mpi_enable}" = "no" -a "${omm_psp_enable}" = "yes"; then
  AC_MSG_ERROR([pspBLAS requires MPI
                  Please adjust your parameters and re-run configure.])
fi

# Clean-up
unset tmp_option
unset tmp_mpi_vars
unset tmp_psp_vars

                    # ------------------------------------ #

#
# Substitute build-system variables
#

# Defaults
AC_SUBST(omm_debug_enable_def)
AC_SUBST(omm_examples_enable_def)
AC_SUBST(omm_randseed_enable_def)
AC_SUBST(omm_mpi_enable_def)
AC_SUBST(omm_mpi_incs_def)
AC_SUBST(omm_mpi_libs_def)
AC_SUBST(omm_linalg_enable_def)
AC_SUBST(omm_linalg_incs_def)
AC_SUBST(omm_linalg_libs_def)
AC_SUBST(omm_psp_enable_def)
AC_SUBST(omm_psp_incs_def)
AC_SUBST(omm_psp_libs_def)
AC_SUBST(omm_msw_enable_def)
AC_SUBST(omm_msw_incs_def)
AC_SUBST(omm_msw_libs_def)

# Triggers
AC_SUBST(omm_debug_enable)
AC_SUBST(omm_examples_enable)
AC_SUBST(omm_randseed_enable)
AC_SUBST(omm_mpi_enable)
AC_SUBST(omm_linalg_enable)
AC_SUBST(omm_psp_enable)
AC_SUBST(omm_msw_enable)

# Initialization types
AC_SUBST(omm_debug_type)
AC_SUBST(omm_examples_type)
AC_SUBST(omm_randseed_type)
AC_SUBST(omm_mpi_type)
AC_SUBST(omm_linalg_type)
AC_SUBST(omm_psp_type)
AC_SUBST(omm_msw_type)

# Internal MPI parameters
AC_SUBST(omm_mpi_cc_set)
AC_SUBST(omm_mpi_cc_wrap)
AC_SUBST(omm_mpi_fc_set)
AC_SUBST(omm_mpi_fc_wrap)

# Build flags for external dependencies
AC_SUBST(omm_mpi_incs)
AC_SUBST(omm_mpi_libs)
AC_SUBST(omm_linalg_incs)
AC_SUBST(omm_linalg_libs)
AC_SUBST(omm_linalg_has_lapack)
AC_SUBST(omm_linalg_has_scalapack)
AC_SUBST(omm_psp_incs)
AC_SUBST(omm_psp_libs)
AC_SUBST(omm_msw_incs)
AC_SUBST(omm_msw_libs)

                    # ------------------------------------ #

#
# Architecture setup
#

# Init MPI parameters, possibly setting and/or wrapping compilers
OMM_MPI_INIT

                    # ------------------------------------ #

#
# C language support
#

# Init C compiler and preprocessor
AC_PROG_CC
AC_PROG_CPP

# C compiler peculiarities (for Libtool)
AM_PROG_CC_C_O

                    # ------------------------------------ #

#
# Fortran language support
#

# Look for the Fortran compiler
if test "${FC}" != "" -a ! -x "${FC}"; then
  omm_fc_probe=`echo "${FC}" | sed -e 's/ .*//'`
  if test ! -x "${omm_fc_probe}"; then
    AC_PATH_PROG([omm_fc_path],[${omm_fc_probe}])
    if test "${omm_fc_path}" = ""; then
      AC_MSG_ERROR([could not run Fortran compiler "${FC}"])
    fi
  fi
fi
AC_PROG_FC

# Fail if no Fortran compiler is available
if test "${FC}" = ""; then
  AC_MSG_ERROR([no Fortran compiler available])
fi

# Set default file extensions
OMM_FC_EXTENSIONS

# Get module file case
OMM_FC_MOD_CASE

                    # ------------------------------------ #

#
# Libtool configuration
#

# Init Libtool (must be done once compilers are fully set)
LT_INIT
LT_PREREQ([2.4.2])
LT_LANG([Fortran])
LTOBJEXT="lo"
AC_SUBST(LTOBJEXT)

                    # ------------------------------------ #

#
# Multicore architecture support
#

# Look for MPI
if test "${omm_mpi_enable}" != "no"; then
  OMM_MPI_DETECT
  if test "${omm_mpi_ok}" = "yes"; then
    AC_DEFINE([MPI], 1, [Define to 1 to enable MPI support.])
    omm_mpi_enable="yes"
  else
    if test "${omm_mpi_enable}" = "auto"; then
      AC_MSG_NOTICE([disabling MPI support])
      omm_mpi_enable="no"
    else
      AC_MSG_FAILURE([MPI support is broken - please check your configuration])
    fi
  fi
fi

# From now on, the MPI trigger must be yes or no
if test "${omm_mpi_enable}" != "no" -a "${omm_mpi_enable}" != "yes"; then
  AC_MSG_ERROR([unexpected MPI trigger value: '${omm_mpi_enable}'])
fi

# Let the test suite know about mpirun
TESTS_CMD="${MPIRUN}"
AC_SUBST(TESTS_CMD)

                    # ------------------------------------ #

#
# External dependencies
#

# Set linear algebra includes
if test "${omm_linalg_type}" != "dir"; then
  if test -n "${LINALG_INCLUDES}"; then
    omm_linalg_incs="${LINALG_INCLUDES}"
  else
    omm_linalg_incs="${omm_linalg_incs_def}"
  fi
fi

# Set linear algebra libraries
if test "${omm_linalg_type}" != "dir"; then
  if test -n "${LINALG_LIBS}"; then
    omm_linalg_libs="${LINALG_LIBS}"
  else
    omm_linalg_libs="${omm_linalg_libs_def}"
  fi
fi

# Look for linear algebra
AC_MSG_CHECKING([whether to look for linear algebra])
AC_MSG_RESULT([${omm_linalg_enable}])
if test "${omm_linalg_enable}" = "yes"; then
  AC_MSG_CHECKING([how linear algebra parameters have been set])
  AC_MSG_RESULT([${omm_linalg_type}])
  AC_MSG_CHECKING([for linear algebra include flags])
  AC_MSG_RESULT([${omm_linalg_incs}])
  AC_MSG_CHECKING([for linear algebra library flags])
  AC_MSG_RESULT([${omm_linalg_libs}])
  OMM_LINALG_DETECT
  if test "${omm_linalg_ok}" != "yes" -a "${omm_linalg_type}" = "def"; then
    AC_MSG_NOTICE([retrying with generic linear algebra library flags (Netlib)])
    LIBS="-lscalapack -lblacs -lblacsCinit -lblacsF77init -llapack -lblas ${LIBS}"
    OMM_LINALG_DETECT
  fi
  if test "${omm_linalg_ok}" = "yes"; then
    CPPFLAGS="${CPPFLAGS} ${omm_linalg_incs}"
    LIBS="${omm_linalg_libs} ${LIBS}"
  else
    AC_MSG_FAILURE([linear algebra is missing or incomplete])
  fi
fi

                              # ---------------- #

# Check that MPI is still enabled
if test "${omm_mpi_enable}" = "no" -a "${omm_psp_enable}" = "yes"; then
  AC_MSG_FAILURE([pspBLAS requires MPI
                  Please adjust your parameters and re-run configure.])
fi

# Set pspBLAS includes
if test "${omm_psp_type}" != "dir"; then
  if test -n "${PSP_INCLUDES}"; then
    omm_psp_incs="${PSP_INCLUDES}"
  else
    omm_psp_incs="${omm_psp_incs_def}"
  fi
fi

# Set pspBLAS libraries
if test "${omm_psp_type}" != "dir"; then
  if test -n "${PSP_LIBS}"; then
    omm_psp_libs="${PSP_LIBS}"
  else
    omm_psp_libs="${omm_psp_libs_def}"
  fi
fi

# Look for pspBLAS
AC_MSG_CHECKING([whether to look for pspBLAS])
AC_MSG_RESULT([${omm_psp_enable}])
if test "${omm_psp_enable}" = "yes"; then
  AC_MSG_CHECKING([how pspBLAS parameters have been set])
  AC_MSG_RESULT([${omm_psp_type}])
  AC_MSG_CHECKING([for pspBLAS include flags])
  AC_MSG_RESULT([${omm_psp_incs}])
  AC_MSG_CHECKING([for pspBLAS library flags])
  AC_MSG_RESULT([${omm_psp_libs}])
  OMM_PSP_DETECT
  if test "${omm_psp_ok}" = "yes"; then
    CPPFLAGS="${CPPFLAGS} ${omm_psp_incs}"
    LIBS="${omm_psp_libs} ${LIBS}"
  else
    AC_MSG_FAILURE([pspBLAS is missing or incomplete])
  fi
fi

# Set MatrixSwitch includes
if test "${omm_msw_type}" != "dir"; then
  if test -n "${MSW_INCLUDES}"; then
    omm_msw_incs="${MSW_INCLUDES}"
  else
    omm_msw_incs="${omm_msw_incs_def}"
  fi
fi

# Set MatrixSwitch libraries
if test "${omm_msw_type}" != "dir"; then
  if test -n "${MSW_LIBS}"; then
    omm_msw_libs="${MSW_LIBS}"
  else
    omm_msw_libs="${omm_msw_libs_def}"
  fi
fi

# Look for MatrixSwitch
AC_MSG_CHECKING([whether to look for MatrixSwitch])
AC_MSG_RESULT([${omm_msw_enable}])
if test "${omm_msw_enable}" = "yes"; then
  AC_MSG_CHECKING([how MatrixSwitch parameters have been set])
  AC_MSG_RESULT([${omm_msw_type}])
  AC_MSG_CHECKING([for MatrixSwitch include flags])
  AC_MSG_RESULT([${omm_msw_incs}])
  AC_MSG_CHECKING([for MatrixSwitch library flags])
  AC_MSG_RESULT([${omm_msw_libs}])
  OMM_MSW_DETECT
  if test "${omm_msw_ok}" = "yes"; then
    CPPFLAGS="${CPPFLAGS} ${omm_msw_incs}"
    LIBS="${omm_msw_libs} ${LIBS}"
  else
    AC_MSG_FAILURE([MatrixSwitch is missing or incomplete])
  fi
fi

                    # ------------------------------------ #

#
# Adjust build parameters according to configure options
#

# Reset option counter
omm_cpp_options=""

# ISO C Bindings are always available nowadays
AC_DEFINE([CBIND], 1, [Define to 1 to use ISO C Bindings])
omm_cpp_options="${omm_cpp_options} CBIND"

# Conditional compilation (alphabetical order)
if test "${omm_debug_enable}" = "yes"; then
  AC_DEFINE([DEBUG_MODE], 1,
    [Define to 1 if you want to enable enhanced debugging information.])
  omm_cpp_options="${omm_cpp_options} DEBUG_MODE"
fi
if test "${omm_randseed_enable}" = "no"; then
  AC_DEFINE([NORAND], 1,
    [Define to 1 if you want a deterministic random seed.])
  omm_cpp_options="${omm_cpp_options} NORAND"
fi
AM_CONDITIONAL([DO_BUILD_EXAMPLES], [test "${enable_examples}" = "yes"])

# External dependencies (alphabetical order)
if test "${omm_linalg_has_lapack}" = "yes"; then
  AC_DEFINE([HAVE_LAPACK], 1, [Define to 1 if you have LAPACK libraries available.])
  omm_cpp_options="${omm_cpp_options} HAVE_LAPACK"
fi
if test "${omm_linalg_has_scalapack}" = "yes"; then
  AC_DEFINE([HAVE_SCALAPACK], 1, [Define to 1 if you have ScaLAPACK libraries available.])
  omm_cpp_options="${omm_cpp_options} HAVE_SCALAPACK"
fi
if test "${omm_psp_ok}" = "yes"; then
  AC_DEFINE([HAVE_PSPBLAS], 1, [Define to 1 if you have the pspBLAS library available.])
  omm_cpp_options="${omm_cpp_options} HAVE_PSPBLAS"
fi

                    # ------------------------------------ #

# Write down YAML configuration now (allows full report with deferred errors)
AC_OUTPUT([config/data/libOMM-config.yaml])

                    # ------------------------------------ #

# Report configuration
AC_MSG_NOTICE([])
AC_MSG_NOTICE([Final build parameters])
AC_MSG_NOTICE([----------------------])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([TSTAMP     = ${omm_timestamp}])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([DEBUG    = ${omm_debug_enable} (init: ${omm_debug_type})])
AC_MSG_NOTICE([EXAMPLES = ${omm_examples_enable} (init: ${omm_examples_type})])
AC_MSG_NOTICE([RANDSEED = ${omm_randseed_enable} (init: ${omm_randseed_type})])
AC_MSG_NOTICE([MPI      = ${omm_mpi_enable} (init: ${omm_mpi_type})])
AC_MSG_NOTICE([LINALG   = ${omm_linalg_enable} (init: ${omm_linalg_type})])
AC_MSG_NOTICE([PSP      = ${omm_psp_enable} (init: ${omm_psp_type})])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([CPP      = ${CPP}])
AC_MSG_NOTICE([CPPFLAGS = ${CPPFLAGS}])
AC_MSG_NOTICE([FC       = ${CC}])
AC_MSG_NOTICE([MPIFC    = ${MPICC}])
AC_MSG_NOTICE([FC       = ${FC}])
AC_MSG_NOTICE([MPIFC    = ${MPIFC}])
AC_MSG_NOTICE([FCFLAGS  = ${FCFLAGS}])
AC_MSG_NOTICE([LDFLAGS  = ${LDFLAGS}])
AC_MSG_NOTICE([LIBS     = ${LIBS}])
AC_MSG_NOTICE([OPTIONS  = ${omm_cpp_options}])
AC_MSG_NOTICE([])

                    # ------------------------------------ #

#
# Output configuration
#

AC_CONFIG_FILES([
  config/data/libOMM.pc
  Makefile
  doc/Makefile
  src/Makefile
  examples/Makefile
])

AC_OUTPUT
