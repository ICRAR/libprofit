#!/bin/bash
#
# Travis CI before-install script
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

# cxxtest can be used from sources, which is the easiest way to go about it
git clone https://github.com/CxxTest/cxxtest
export CMAKE_PREFIX_PATH=$PWD/cxxtest

# In MacOS we simply need to brew install some things
if [ "${TRAVIS_OS_NAME}" = "osx" ]
then

	# The xcode7.3 osx image needs an update
	if [ "${XCODE}" = "7.3" ]
	then
		brew update
	fi

	# Minimal dependencies for testing
	brew install gsl

	# 7.3 (El Capitan) mostly doesn't have associated bottles anymore,
	# and fftw is pulling both gcc and open-mpi, which take way too long
	# to compile; let's simply skip this one
	if [ "${XCODE}" != "7.3" ]; then
		brew install fftw
	fi
	return
fi

# Otherwise this is a linux box
# We use cmake, but in trusty we still have 2.8...
curl -O https://cmake.org/files/v3.1/cmake-3.1.3-Linux-x86_64.tar.gz
tar -xf cmake-3.1.3-Linux-x86_64.tar.gz
export PATH=${TRAVIS_BUILD_DIR}/cmake-3.1.3-Linux-x86_64/bin:$PATH
export CMAKE_MODULE_PATH=${TRAVIS_BUILD_DIR}/cmake-3.1.3-Linux-x86_64/share/cmake-3.1/Modules:${CMAKE_MODULE_PATH}

# We are collecting coverage
if [ "$COMPILER" == "g++-6" ];
then
	pip install --user cpp-coveralls
fi
