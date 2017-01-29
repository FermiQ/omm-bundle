check_PROGRAMS = \
  example_gamma.x \
  example_kpoints.x \
  example_gamma_wrapper.x \
  example_kpoints_wrapper.x

example_gamma_x_SOURCES = example_gamma.F90
example_gamma_x_FCFLAGS = -I../src -DSRCDIR=\'$(srcdir)\'
example_gamma_x_LDADD = ../src/libomm.la

example_kpoints_x_SOURCES = example_kpoints.F90
example_kpoints_x_FCFLAGS = -I../src -DSRCDIR=\'$(srcdir)\'
example_kpoints_x_LDADD = ../src/libomm.la

example_gamma_wrapper_x_SOURCES = example_gamma_wrapper.F90
example_gamma_wrapper_x_FCFLAGS = -I../src -DSRCDIR=\'$(srcdir)\'
example_gamma_wrapper_x_LDADD = ../src/libomm.la

example_kpoints_wrapper_x_SOURCES = example_kpoints_wrapper.F90
example_kpoints_wrapper_x_FCFLAGS = -I../src -DSRCDIR=\'$(srcdir)\'
example_kpoints_wrapper_x_LDADD = ../src/libomm.la

LOG_COMPILER = $(TESTS_CMD)

TESTS = \
  example_gamma.x \
  example_kpoints.x \
  example_gamma_wrapper.x \
  example_kpoints_wrapper.x

CLEANFILES = \
  libOMM.log
