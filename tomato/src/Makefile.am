# -*- Automake -*-
#
# Makefile for the Tomato package
#
# Copyright (C) 2016 Yann Pouillon
#
# This file is part of the Tomato software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

                    # ------------------------------------ #

#
# Main source files
#

# Common source files
tom_srcs = \
  tomato.F90

# Fortran modules
# TODO: support both lower-case and upper-case file names
tom_f03_mods = \
  tomato.$(MODEXT)

# Libraries to install
lib_LTLIBRARIES = libtomato.la

libtomato_la_SOURCES = $(tom_srcs)
libtomato_la_LDFLAGS = -version-info 0:0:0

# Fortran modules to install
f03moddir = $(includedir)
install-data-local:
	$(INSTALL) -d -m 755 $(DESTDIR)$(f03moddir)
	$(INSTALL) -m 644 $(tom_f03_mods) $(DESTDIR)$(f03moddir)

uninstall-local:
	cd $(DESTDIR)$(f03moddir) && rm -f $(tom_f03_mods)

# Local cleaning
clean-local:
	rm -f *.mod *.MOD *.obj

# Files to mark explicitly for distribution
EXTRA_DIST = tomato_TB.F90
