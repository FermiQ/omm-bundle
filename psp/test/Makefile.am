check_PROGRAMS = pdgemmScaling.x pzgemmScaling.x

pdgemmScaling_x_SOURCES = pdgemmScaling.F90
pdgemmScaling_x_FCFLAGS = -I../src
pdgemmScaling_x_LDADD = ../src/libpspblas.la

pzgemmScaling_x_SOURCES = pzgemmScaling.F90
pzgemmScaling_x_FCFLAGS = -I../src
pzgemmScaling_x_LDADD = ../src/libpspblas.la

LOG_COMPILER = $(TESTS_CMD)

TESTS = \
  pdgemmScaling.x \
  pzgemmScaling.x

