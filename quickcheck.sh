#!/bin/sh
#
# Copyright (C) 2017 Yann Pouillon <notifications@materialsevolution.es>
#
# This file is part of libOMM.
#
# libOMM is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, version 3 of the License, or (at your option) any later
# version.
#
# libOMM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with libOMM.  If not, see <http://www.gnu.org/licenses/> or write to the
# Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
# 02110-1301  USA.

# Note: this script is temporary and will be removed upon release.

# Stop at first error and echo commands
set -ev

# Check that we are in the correct directory
test -x "pspBLAS/quickcheck.sh" -a \
     -x "MatrixSwitch/quickcheck.sh" -a \
     -x "libOMM/quickcheck.sh" || exit 0

# Prepare source tree
./wipeout.sh
./autogen.sh

# Run individual scripts
for srcdir in pspBLAS MatrixSwitch libOMM; do
  cd "${srcdir}"
  ./quickcheck.sh
  cd ..
done
