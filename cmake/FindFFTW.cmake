#
# FFTW library finding routine for cmake
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2017
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
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

find_path(FFTW_INCLUDE_DIR fftw3.h)
find_library(FFTW_LIBRARIES NAMES fftw3)

# If OpenMP has been found already we try to include the OpenMP flavour of FFTW
if (OPENMP_FOUND)
	find_library(FFTW_OPENMP_LIB NAMES fftw3_omp)
	if (EXISTS ${FFTW_OPENMP_LIB})
		set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_OPENMP_LIB})
		set(FFTW_OPENMP_FOUND ON)
	endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDE_DIR)
