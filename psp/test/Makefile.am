check_PROGRAMS = pdgemmScaling.x \
  pzgemmScaling.x \
  testMPI.x \
  testMPI2.x \
  testMPI3.x \
  testList.x \
  spBLAS.x \
  pdgemm.x \
  matSum.x

pdgemmScaling_x_SOURCES = pdgemmScaling.F90
pdgemmScaling_x_FCFLAGS = -I../src
pdgemmScaling_x_LDADD = ../src/libpspblas.la

pzgemmScaling_x_SOURCES = pzgemmScaling.F90
pzgemmScaling_x_FCFLAGS = -I../src
pzgemmScaling_x_LDADD = ../src/libpspblas.la

testMPI_x_SOURCES = testMPI.F90
testMPI_x_FCFLAGS = -I../src
testMPI_x_LDADD = ../src/libpspblas.la

testMPI2_x_SOURCES = testMPI2.F90
testMPI2_x_FCFLAGS = -I../src
testMPI2_x_LDADD = ../src/libpspblas.la

testMPI3_x_SOURCES = testMPI3.F90
testMPI3_x_FCFLAGS = -I../src
testMPI3_x_LDADD = ../src/libpspblas.la

testList_x_SOURCES = testList.F90
testList_x_FCFLAGS = -I../src
testList_x_LDADD = ../src/libpspblas.la

spBLAS_x_SOURCES = spBLAS.F90
spBLAS_x_FCFLAGS = -I../src
spBLAS_x_LDADD = ../src/libpspblas.la

pdgemm_x_SOURCES = pdgemm.F90
pdgemm_x_FCFLAGS = -I../src
pdgemm_x_LDADD = ../src/libpspblas.la

matSum_x_SOURCES = matSum.F90
matSum_x_FCFLAGS = -I../src
matSum_x_LDADD = ../src/libpspblas.la

LOG_COMPILER = $(TESTS_CMD)

TESTS = \
  pdgemmScaling.x \
  pzgemmScaling.x \
  testMPI.x \
  testMPI2.x \
  testMPI3.x \
  testList.x \
  spBLAS.x \
  pdgemm.x \
  matSum.x

CLEANFILES = \
  zgemm.txt \
  dgemm.txt \
  zlap.txt \
  zspmm.txt \
  dspmm.txt \
  zmspm.txt \
  dlap.txt \
  dmspm.txt

