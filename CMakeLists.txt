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
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMakeModules/)

# Macros to check for presence of GSL and R/Rmath
macro(find_gsl)
	find_package(GSL)
	if( GSL_FOUND )
		add_definitions(-DHAVE_GSL)
		include_directories(${GSL_INCLUDE_DIRS})
		set(profit_LIBS ${GSL_LIBRARIES})
	endif()
endmacro()

macro(find_r_and_rmath)
	find_package(PkgConfig)
	if( PKG_CONFIG_FOUND )
		pkg_check_modules(R libR)
		pkg_check_modules(RMATH libRmath)
		if( R_FOUND AND RMATH_FOUND )
			add_definitions(-DHAVE_R -DMATHLIB_STANDALONE)
			include_directories(${R_INCLUDE_DIRS} ${RMATH_INCLUDE_DIRS})
			set(profit_LIBS ${R_LIBRARIES} ${RMATH_LIBRARIES})
		endif( R_FOUND AND RMATH_FOUND )
	endif( PKG_CONFIG_FOUND )
endmacro()

# Depending on the user's preference we check GSL or R/Rmath first
# We default to use GSL
if( LIBPROFIT_USE_R )
	find_r_and_rmath()
	if( NOT R_FOUND OR NOT RMATH_FOUND)
		find_gsl()
	endif()
else()
	find_gsl()
	if( NOT GSL_FOUND )
		find_r_and_rmath()
	endif()
endif()


# Common definitions for the shared lib and the binary
add_definitions(-DPROFIT_VERSION="${PROFIT_VERSION}")
include_directories(${CMAKE_SOURCE_DIR})

# The shared library
set(PROFIT_SRC
   src/convolve.cpp
   src/coresersic.cpp
   src/ferrer.cpp
   src/king.cpp
   src/moffat.cpp
   src/profit.cpp
   src/psf.cpp
   src/radial.cpp
   src/sersic.cpp
   src/sky.cpp
   src/utils.cpp
)
add_library(profit SHARED ${PROFIT_SRC})
target_link_libraries(profit ${profit_LIBS})

# The executable
add_executable(profit-cli src/profit-cli.cpp)
target_link_libraries(profit-cli profit)

# Installing stuff
install(TARGETS profit profit-cli
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/profit DESTINATION include)


# Tests
if( LIBPROFIT_TEST )
	find_package(CxxTest)
	if( CXXTEST_FOUND )
		include_directories(${CXXTEST_INCLUDE_DIR})
		enable_testing()
		add_subdirectory(tests)
	endif()
endif()
