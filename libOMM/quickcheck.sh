#!/bin/sh
#
# Copyright (C) 2016 Yann Pouillon <notifications@materialsevolution.es>
#
# This file is part of LibOMM.
#
# LibOMM is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, version 3 of the License, or (at your option) any later
# version.
#
# LibOMM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with LibOMM.  If not, see <http://www.gnu.org/licenses/> or write to the
# Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
# 02110-1301  USA.

# Note: this script is temporary and will be removed upon release.

# Stop at first error and echo commands
set -ev

# Check that we are in the correct directory
test -s "configure.ac" -a -s "src/omm.F90" || exit 0

# Required install dirs
MSW_ROOT="${PWD}/../tmp-msw"
PSP_ROOT="${PWD}/../tmp-psp"

# Init build parameters
export CC="mpicc"
export FC="mpif90"
export CFLAGS="-O0 -g3 -ggdb -Wall -Wextra -fbounds-check -fno-inline"
export FCFLAGS="-O0 -g3 -ggdb -Wall -Wextra -fbounds-check -fno-inline"
export MSW_INCLUDES="-I${MSW_ROOT}/include"
export MSW_LIBS="-L${MSW_ROOT}/lib -lMatrixSwitch"
export PSP_INCLUDES="-I${PSP_ROOT}/include"
export PSP_LIBS="-L${PSP_ROOT}/lib -lpspblas"

# Prepare source tree
./wipeout.sh
./autogen.sh

# Check default build
mkdir tmp-minimal
cd tmp-minimal
../configure
make dist
make
make clean && make -j4
make check
mkdir install-minimal
make install DESTDIR="${PWD}/install-minimal"
ls -lR install-minimal >install-minimal.log
cd ..

# Check examples build
mkdir tmp-examples
cd tmp-examples
../configure --enable-examples
make -j4
make check -j4
cd ..

# Check Linalg build
mkdir tmp-linalg
cd tmp-linalg
../configure \
  LINALG_LIBS="-llapack -lblas"
make -j4
make check -j4
cd ..

# Check bare MPI build
mkdir tmp-mpi
cd tmp-mpi
../configure --enable-examples CC="mpicc" FC="mpif90"
make -j4
make check -j4
cd ..

# Check MPI + Linalg build
mkdir tmp-mpi-linalg
cd tmp-mpi-linalg
../configure \
  LINALG_LIBS="-lscalapack -lblacs -lblacsCinit -lblacsF77init -llapack -lblas" \
  --enable-examples \
  CC="mpicc" FC="mpif90"
make -j4
make check -j4
cd ..

# Check MPI + Psp build
#mkdir tmp-mpi-psp
#cd tmp-mpi-psp
#../configure \
#  --with-psp="${PWD}/../../tmp-psp" \
#  LINALG_LIBS="-lscalapack -lblacs -lblacsCinit -lblacsF77init -llapack -lblas" \
#  --enable-examples \
#  CC="mpicc" FC="mpif90"
#make -j4
#make check -j4
#cd ..

# Make distcheck
mkdir tmp-distcheck
cd tmp-distcheck
../configure
make distcheck -j4
make distcleancheck

# Clean-up the mess
cd ..
rm -rf tmp-minimal tmp-examples tmp-linalg tmp-mpi tmp-mpi-linalg tmp-mpi-psp tmp-distcheck
