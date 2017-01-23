check_PROGRAMS = \
  example_gamma_unpolarized.x \
  example_gamma_unpolarized_overlap.x \
  example_gamma_unpolarized_wrapper.x

example_gamma_unpolarized_x_SOURCES = example_gamma_unpolarized.F90
example_gamma_unpolarized_x_FCFLAGS = -I../src
example_gamma_unpolarized_x_LDADD = ../src/libomm.la

example_gamma_unpolarized_overlap_x_SOURCES = example_gamma_unpolarized_overlap.F90
example_gamma_unpolarized_overlap_x_FCFLAGS = -I../src
example_gamma_unpolarized_overlap_x_LDADD = ../src/libomm.la

example_gamma_unpolarized_wrapper_x_SOURCES = example_gamma_unpolarized.F90
example_gamma_unpolarized_wrapper_x_FCFLAGS = -I../src
example_gamma_unpolarized_wrapper_x_LDADD = ../src/libomm.la

LOG_COMPILER = $(TESTS_CMD)

TESTS = \
  example_gamma_unpolarized.x \
  example_gamma_unpolarized_overlap.x \
  example_gamma_unpolarized_wrapper.x
