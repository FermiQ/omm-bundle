#
# Makefile for the documentation of Libpsp
#

# Markdown sources
psp_mddocs = \
  configuring-libpsp.md \
  contributing-to-libpsp.md \
  hacking-the-build-system.md

# HTML documents
psp_htmldocs = \
  configuring-libpsp.html \
  contributing-to-libpsp.html \
  hacking-the-build-system.html

                    ########################################

# User-defined rules
MARKDOWN = $(top_srcdir)/config/scripts/markdown
MDFLAGS  = 

SUFFIXES = .md .html
.md.html:
	sed -e "s/@TITLE@/$$(head -n 1 $<)/" \
	  $(srcdir)/data/psp_header.html >$@
	$(MARKDOWN) $(MDFLAGS) $< >>$@
	cat $(srcdir)/data/psp_footer.html >>$@

# Targets
htmldocdir = $(docdir)/html
htmldoc_DATA = $(psp_htmldocs)

# Inform Automake
CLEANFILES = $(psp_htmldocs)
EXTRA_DIST = data/psp_header.html data/psp_footer.html $(psp_mddocs)
