# CMakeLists.txt
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2016
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of libprofit.
#
# libprofit is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# libprofit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with libprofit.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.1)

# We have the version here
file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION PROFIT_VERSION)

# Our project...
project(libprofit VERSION ${PROFIT_VERSION} LANGUAGES CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMakeModules/)

# Check if there's GSL or R support, in that order
find_package(GSL)
if( GSL_FOUND )
	add_definitions(-DHAVE_GSL)
	include_directories(${GSL_INCLUDE_DIRS})
	set(profit_LIBS ${GSL_LIBRARIES})
else ()
	find_package(LibR)
	if( LIBR_FOUND )
		add_definitions(-DHAVE_R)
		include_directories(${LIBR_INCLUDE_DIRS})
		set(profit_LIBS ${LIBR_LIBRARIES})
	endif( LIBR_FOUND )
endif( GSL_FOUND )

# Common definitions for the shared lib and the binary
add_definitions(-DPROFIT_VERSION="${PROFIT_VERSION}")
include_directories(include)

# The shared library
set(PROFIT_SRC
   src/analytic.cpp
   src/convolve.cpp
   src/ferrer.cpp
   src/moffat.cpp
   src/profit.cpp
   src/psf.cpp
   src/sersic.cpp
   src/sky.cpp
   src/utils.cpp
)
add_library(profit SHARED ${PROFIT_SRC})
target_link_libraries(profit ${profit_LIBS})
target_compile_features(profit PUBLIC cxx_auto_type cxx_range_for)

# The executable
add_executable(profit-cli src/profit-cli.cpp)
target_link_libraries(profit-cli profit)

# Installing stuff
install(TARGETS profit profit-cli
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)
