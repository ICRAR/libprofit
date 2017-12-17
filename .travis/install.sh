#!/bin/bash
#
# Travis CI install script
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2017
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307  USA
#

cd ${TRAVIS_BUILD_DIR}
mkdir build
cd build

LIBPROFIT_CMAKE_OPTIONS="-DCMAKE_CXX_COMPILER=$COMPILER -DLIBPROFIT_TEST=ON"

# Poor builds disable everything
if [ "${POOR_BUILD}" = "yes" ]
then
	LIBPROFIT_CMAKE_OPTIONS="$LIBPROFIT_CMAKE_OPTIONS -DLIBPROFIT_NO_OPENCL=ON -DLIBPROFIT_NO_OPENMP=ON -DLIBPROFIT_NO_FFTW=ON"
fi

# Go, go, go!
cmake .. ${LIBPROFIT_CMAKE_OPTIONS}
make all
