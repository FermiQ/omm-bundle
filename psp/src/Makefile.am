# -*- Automake -*-
#
# Makefile for the pspBLAS package
#
# Copyright (C) 2016 Yann Pouillon
#
# This file is part of the pspBLAS software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

                    # ------------------------------------ #

#
# Main source files
#

# Common source files
psp_srcs = \
  pspBLAS.F90 \
  pspLevel1/pspLevel1.F90 \
  pspLevel2/pspLevel2.F90 \
  pspLevel3/pspGemm.F90 \
  pspLevel3/pspLevel3.F90 \
  pspLevel3/pspMspm.F90 \
  pspLevel3/pspSpmm.F90 \
  pspMPI/pspMPI.F90 \
  pspUtility/pspUtility.F90 \
  pspVariable/pspVariable.F90

# Fortran modules
# TODO: support both lower-case and upper-case file names
psp_f03_mods = \
  psplevel3.$(MODEXT) \
  pspgemm.$(MODEXT) \
  pspspmm.$(MODEXT) \
  pspmspm.$(MODEXT) \
  pspmpi.$(MODEXT) \
  psplevel2.$(MODEXT) \
  pspblas.$(MODEXT) \
  psputility.$(MODEXT) \
  pspvariable.$(MODEXT) \
  psplevel1.$(MODEXT)

# Libraries to install
lib_LTLIBRARIES = libpspblas.la

libpspblas_la_SOURCES = $(psp_srcs)
libpspblas_la_LDFLAGS = -version-info 0:0:0

# Fortran modules to install
f03moddir = $(includedir)
install-data-local:
	$(INSTALL) -d -m 755 $(DESTDIR)$(f03moddir)
	$(INSTALL) -m 644 $(psp_f03_mods) $(DESTDIR)$(f03moddir)

uninstall-local:
	cd $(DESTDIR)$(f03moddir) && rm -f $(psp_f03_mods)

# Local cleaning
clean-local:
	rm -f *.mod *.MOD *.obj

# Explicit dependencies
include $(srcdir)/dependencies.mf

# Files to mark explicitly for distribution
EXTRA_DIST = dependencies.mf
