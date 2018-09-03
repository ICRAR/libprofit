#!/bin/bash
#
# Travis CI after-success script
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

# We only execute further actions if running a coverage build
if [ "$COMPILER" != "g++-6" ]
then
	return
fi

cd ${TRAVIS_BUILD_DIR}

# Get coverage results and upload to coveralls
coveralls --gcov `which gcov-6` -b build --gcov-options '\-lp' -r . -i ./src -i ./profit -e src/profit-cli.cpp -E '.*fits_utils.*'

# Re-generate the coverage results (manually this time),
# check the code quality, and upload all results to sonarqube
mkdir gcov-reports
cd gcov-reports
gcov-6 -lp ../build/CMakeFiles/profit.dir/src/*.gcno ../build/CMakeFiles/profit-cli.dir/src/*.gcno
cd ..
sonar-scanner \
    -Dsonar.sources=src,profit \
    -Dsonar.projectKey=libprofit \
    -Dsonar.organization=rtobar-github \
    -Dsonar.cfamily.build-wrapper-output=bw-output \
    -Dsonar.exclusions=profit/cl/cl2.hpp \
    -Dsonar.cfamily.gcov.reportsPath=gcov-reports
