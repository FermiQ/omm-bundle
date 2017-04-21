#!/bin/sh
#
# Copyright (C) 2016 Yann Pouillon <notifications@materialsevolution.es>
#
# This file is part of MatrixSwitch.
#
# MatrixSwitch is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, version 3 of the License, or (at your option) any later
# version.
#
# MatrixSwitch is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MatrixSwitch.  If not, see <http://www.gnu.org/licenses/> or write
# to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
# 02110-1301  USA.

# Note: this script is temporary and will be removed upon release.

# Stop at first error and echo commands
set -ev

# Check that we are in the correct directory
test -s "configure.ac" -a -s "src/MatrixSwitch.F90" || exit 0

# Set number of processors for parallel builds (make -j)
make_nprocs="8"

# Init build parameters
export CFLAGS="-O0 -g3 -ggdb -Wall -Wextra -fbounds-check -fno-inline"
export FCFLAGS="-O0 -g3 -ggdb -Wall -Wextra -fbounds-check -fno-inline"

# Prepare source tree
./wipeout.sh
./autogen.sh

# Check default build
mkdir tmp-minimal
cd tmp-minimal
instdir="${PWD}/install-minimal"
../configure \
  --prefix="${instdir}" \
  --disable-debug \
  --without-mpi
sleep 3
make dist
make
make clean && make -j${make_nprocs}
make -j${make_nprocs} check
make -j${make_nprocs} install
ls -lR "${instdir}" >../install-minimal.tmp
cat ../install-minimal.tmp
sleep 3
cd ..

# Check default build
mkdir tmp-default
cd tmp-default
../configure
sleep 3
make -j${make_nprocs}
make -j${make_nprocs} check
sleep 3
cd ..

# Check Linalg build (EasyBuild)
mkdir tmp-linalg1
cd tmp-linalg1
../configure \
  --with-linalg
sleep 3
make -j${make_nprocs}
make -j${make_nprocs} check
sleep 3
cd ..

# Check Linalg build (system libs)
mkdir tmp-linalg2
cd tmp-linalg2
../configure \
  LINALG_LIBS="-llapack -lblas"
sleep 3
make -j${make_nprocs}
make -j${make_nprocs} check
sleep 3
cd ..

# Check bare MPI build
mkdir tmp-mpi
cd tmp-mpi
../configure \
  CC="mpicc" \
  FC="mpif90"
sleep 3
make -j${make_nprocs}
make -j${make_nprocs} check
sleep 3
cd ..

# Check MPI + Linalg build
mkdir tmp-mpi-linalg
cd tmp-mpi-linalg
../configure \
  LINALG_LIBS="-lscalapack -lopenblas" \
  CC="mpicc" \
  FC="mpifort"
sleep 3
make -j${make_nprocs}
make -j${make_nprocs} check
sleep 3
cd ..

# Check MPI + pspBLAS build
mkdir tmp-mpi-psp
cd tmp-mpi-psp
if test -s "../../build-omm"; then
  instdir="${PWD}/../../tmp-matrixswitch"
else
  instdir="${PWD}/tmp-install-all"
fi
../configure \
  --prefix="${instdir}" \
  --with-psp="${PWD}/../../tmp-pspblas" \
  LINALG_LIBS="-lscalapack -lopenblas" \
  CC="mpicc" \
  FC="mpifort"
sleep 3
make -j${make_nprocs}
make -j${make_nprocs} check
sleep 3
make -j${make_nprocs} install
cd ..

# Make distcheck
mkdir tmp-distcheck
cd tmp-distcheck
../configure
sleep 3
make -j${make_nprocs} distcheck
make distcleancheck

# Clean-up the mess
cd ..
rm -rf tmp-minimal tmp-default tmp-linalg1 tmp-linalg2 tmp-mpi tmp-mpi-linalg tmp-mpi-psp tmp-distcheck
