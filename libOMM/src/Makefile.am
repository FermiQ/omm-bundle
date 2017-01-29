# -*- Automake -*-
#
# Makefile for the libOMM package
#
# Copyright (C) 2016 Yann Pouillon
#
# This file is part of the libOMM software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

                    # ------------------------------------ #

#
# Main source files
#

# Common source files
omm_srcs = \
  omm.F90 \
  omm_callback.F90 \
  omm_ops.F90 \
  omm_params.F90 \
  omm_quartic.F90 \
  omm_rand.F90 \
  omm_wrapper.F90

# Fortran modules
# TODO: support both lower-case and upper-case file names
omm_f03_mods = \
  omm_ops.$(MODEXT) \
  omm_params.$(MODEXT) \
  omm_rand.$(MODEXT)

# Libraries to install
lib_LTLIBRARIES = libOMM.la

libOMM_la_SOURCES = $(omm_srcs)
libOMM_la_LDFLAGS = -version-info 0:0:0

# Fortran modules to install
f03moddir = $(includedir)
install-data-local:
	$(INSTALL) -d -m 755 $(DESTDIR)$(f03moddir)
	$(INSTALL) -m 644 $(omm_f03_mods) $(DESTDIR)$(f03moddir)

uninstall-local:
	cd $(DESTDIR)$(f03moddir) && rm -f $(omm_f03_mods)

# Local cleaning
clean-local:
	rm -f *.mod *.MOD *.obj

# Explicit dependencies
include $(srcdir)/dependencies.mf

# Files to mark explicitly for distribution
EXTRA_DIST = dependencies.mf
