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

# Extract version information from VERSION file
file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION PROFIT_VERSION)
message("-- Building libproft version ${PROFIT_VERSION}")
string(REPLACE "." ";" VERSION_LIST ${PROFIT_VERSION})
string(REPLACE "-" ";" VERSION_LIST "${VERSION_LIST}")

list(LENGTH VERSION_LIST PROFIT_VERSION_N_ITEMS)
list(GET VERSION_LIST 0 PROFIT_VERSION_MAJOR)
list(GET VERSION_LIST 1 PROFIT_VERSION_MINOR)
list(GET VERSION_LIST 2 PROFIT_VERSION_PATCH)
if (${PROFIT_VERSION_N_ITEMS} GREATER 3)
	list(GET VERSION_LIST 3 PROFIT_VERSION_SUFFIX)
else ()
	set(PROFIT_VERSION_SUFFIX "")
endif()

# We have only C++ sources, but until cmake 3.4 the check_symbol_exists
# macro didn't work properly if only defining C++ as the project language
if( ${CMAKE_MAJOR_VERSION} EQUAL 3 AND ${CMAKE_MINOR_VERSION} LESS 4 )
	set(PROFIT_LANGS C CXX)
else()
	set(PROFIT_LANGS CXX)
endif()

# Our project...
project(libprofit VERSION "${PROFIT_VERSION_MAJOR}.${PROFIT_VERSION_MINOR}.${PROFIT_VERSION_PATCH}" LANGUAGES ${PROFIT_LANGS})
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/)

# Compile in release mode by default
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif()

# Options users can give on the command line
option(LIBPROFIT_USE_R     "Try to compile libprofit first against R/Rmath, and then with GSL" OFF)
option(LIBPROFIT_TEST      "Include test compilation in the build" OFF)
option(LIBPROFIT_DEBUG     "Include some debugging code in libprofit" OFF)
option(LIBPROFIT_NO_OPENCL "Don't attempt to include OpenCL support in libprofit" OFF)
option(LIBPROFIT_NO_OPENMP "Don't attempt to include OpenMP support in libprofit" OFF)
option(LIBPROFIT_NO_FFTW   "Don't attempt to include FFTW support in libprofit" OFF)
option(LIBPROFIT_NO_SIMD   "Don't attempt to use SIMD (SSE2/AVX) extensions in libprofit" OFF)

#
# Macros to check for the presence of:
#  * GSL
#  * R/Rmath
#  * OpenCL
#  * OpenMP
#  * FFTW
#
# From those, only GSL or R/Rmath are required (one of the two), the rest is
# optional.
#
macro(find_gsl)

	# Try with FindGSL.cmake first, if found
	# Otherwise try with pkg-config
	find_package(GSL)
	if( NOT GSL_FOUND )
		find_package(PkgConfig)
		if( PKG_CONFIG_FOUND )
			pkg_check_modules(GSL gsl)
		endif()
	endif()

	if( GSL_FOUND )
		set(PROFIT_USES_GSL ON)
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
			set(PROFIT_USES_R ON)
			add_definitions(-DMATHLIB_STANDALONE)
			include_directories(${R_INCLUDE_DIRS} ${RMATH_INCLUDE_DIRS})
			set(profit_LIBS ${R_LIBRARIES} ${RMATH_LIBRARIES})
			link_directories(${R_LIBRARY_DIRS} ${RMATH_LIBRARY_DIRS})
		endif( R_FOUND AND RMATH_FOUND )
	endif( PKG_CONFIG_FOUND )
endmacro()

macro(find_opencl)
	find_package(OpenCL)
	if( OpenCL_FOUND )
		if ("${OpenCL_VERSION_STRING}" STREQUAL "")
			message("-- OpenCL found but no version reported. Compiling without OpenCL support")
		else ()
			include_directories(${OpenCL_INCLUDE_DIRS})
			set(PROFIT_OPENCL ON)
			set(PROFIT_OPENCL_MAJOR ${OpenCL_VERSION_MAJOR})
			set(PROFIT_OPENCL_MINOR ${OpenCL_VERSION_MINOR})

			# We usually set our target version to be OpenCL 1.1,
			# but in some platforms we know we can aim higher
			# (and save use lots of warnings)
			set(PROFIT_OPENCL_TARGET_VERSION 110)
			if (APPLE)
				string(REGEX REPLACE "^([0-9]+)\\.([0-9]+).*$" "\\1"  DARWIN_VERSION_MAJOR "${CMAKE_SYSTEM_VERSION}")
				if ("${DARWIN_VERSION_MAJOR}" GREATER 11)
					set(PROFIT_OPENCL_TARGET_VERSION 120)
				endif()
			endif()

			# Some platforms (Ubuntu 18.04, for instance) are already
			# shipping with OpenCL versions above 2.0, which is the
			# maximum version our bundled cl2.hpp header supports
			# (and there doesn't seem to be many updates upstream).
			# In such cases, the header will complain about the
			# version not having a supported value, and will default
			# to using 2.0. This is fine (we use up to 1.2 features),
			# but it's nicer to *not* get a warning because of that.
			if (${PROFIT_OPENCL_MAJOR} EQUAL 2 AND NOT ${PROFIT_OPENCL_MINOR} EQUAL 0)
				set(PROFIT_OPENCL_MINOR 0)
			endif()
			set(profit_LIBS ${profit_LIBS} ${OpenCL_LIBRARIES})
		endif()
	endif()
endmacro()

macro(find_openmp)
	find_package(OpenMP)
	if( OPENMP_FOUND )

		# We require at least OpenMP 3.0, so let's check for that
		set(OPENMP_VERSION_CHECK_SOURCE "
#include <stdio.h>

int main(int argc, char *argv[]) {
#if _OPENMP >= 200203
	return 0;
#else
	fail to compile please
#endif
}")
		set(WORK_DIR ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/OpenMPVersionCheck)
		set(SRC_FILE ${WORK_DIR}/ompver.cpp)
		file(WRITE ${SRC_FILE} "${OPENMP_VERSION_CHECK_SOURCE}")
		try_compile(COMPILE_RESULT ${CMAKE_BINARY_DIR} ${SRC_FILE}
		            COMPILE_DEFINITIONS ${OpenMP_CXX_FLAGS})

		if (COMPILE_RESULT)
			set(PROFIT_OPENMP ON)
		else()
			message("-- OpenMP found, but is not >= 2.0. Compiling without OpenMP support")
		endif()
	endif()

	if (PROFIT_OPENMP)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()
endmacro()

macro(find_fftw)
	find_package(FFTW)
	if( FFTW_FOUND )
		set(PROFIT_FFTW ON)
		if( FFTW_OPENMP_FOUND )
			set(PROFIT_FFTW_OPENMP ON)
		endif()
		include_directories(${FFTW_INCLUDE_DIR})
		set(profit_LIBS ${profit_LIBS} ${FFTW_LIBRARIES})
	endif()
endmacro(find_fftw)

macro(find_simd_extensions)
	set(_SSE2 "SSE2;emmintrin.h;__m128d;_mm_setzero_pd;_mm_add_pd")
	set(_AVX "AVX;immintrin.h;__m256d;_mm256_setzero_pd;_mm256_add_pd")
	set(WORK_DIR ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/IntrinsicsCheck)
	foreach(_instruction_details "${_SSE2}" "${_AVX}")
		list(GET _instruction_details 0 _name)
		list(GET _instruction_details 1 _include)
		list(GET _instruction_details 2 _type)
		list(GET _instruction_details 3 _zero)
		list(GET _instruction_details 4 _add)

		set(_src "
#include <${_include}>
int main(int argc, char *argv[]) {
	volatile ${_type} x = ${_zero}();
	volatile ${_type} y = ${_add}(x, x);
}")
		set(SRC_FILE ${WORK_DIR}/${_name}.cpp)
		file(WRITE ${SRC_FILE} "${_src}")
		try_compile(PROFIT_HAS_${_name} ${CMAKE_BINARY_DIR} ${SRC_FILE})

		if (PROFIT_HAS_${_name})
			message(STATUS "Found ${_name} support")
		else()
			message(STATUS "No ${_name} support found")
		endif()
	endforeach()
endmacro()

#
# Depending on the user's preference we check GSL or R/Rmath first
# We default to use GSL
#
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

#
# None found, what an ill luck
if( NOT GSL_FOUND )
	if( NOT R_FOUND OR NOT RMATH_FOUND )
		message(FATAL_ERROR "Neither GSL nor R/Rmath were found, cannot continue")
	endif()
endif()

#
# Check if there's OpenCL (users might opt out)
#
if( NOT LIBPROFIT_NO_OPENCL )
	find_opencl()
endif()

#
# Check if there's OpenMP support (users might opt out)
#
if( NOT LIBPROFIT_NO_OPENMP )
	find_openmp()
endif()

#
# Check if there's FFTW support (users might opt out)
#
if( NOT LIBPROFIT_NO_FFTW )
	find_fftw()
endif()

#
# Include debugging stuff? (off by default)
#
if( LIBPROFIT_DEBUG )
	set(PROFIT_DEBUG ON)
endif()

#
# Check whether we can use CPU extensions
#
if (NOT LIBPROFIT_NO_SIMD)
	find_simd_extensions()
endif()

#
# Create the final configuration file that records
# the compile-time options used for this libprofit build
#
configure_file(include/profit/config.h.in profit/config.h @ONLY)

#
# Common definitions for the shared lib and the binary
#
include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_BINARY_DIR})

#
# MSVC builds need the 3rdparty/getopt submodule to be checked out,
# and added to the list of include directories (so getopt.h is found)
#
if (MSVC)
	set(PROFIT_GETOPT_DIR "${CMAKE_SOURCE_DIR}/3rdparty/getopt")
	if (NOT EXISTS "${PROFIT_GETOPT_DIR}")
		message(FATAL_ERROR
		        "git submodule '3rdparty/getopt' is not checked out but is needed for Windows builds")
	endif()
	include_directories("${PROFIT_GETOPT_DIR}")
endif()

# Compiler-specific tweaks
if (MSVC)
	# Do not warn about unsecure standard library calls
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

#
# Generate header files with the OpenCL kernel sources out of each individual
# .cl file
#
if (PROFIT_OPENCL)
set(OPENCL_KERNEL_HEADERS "")
	foreach(_cl_kernel brokenexponential common convolve coresersic ferrer king moffat sersic)
		foreach(_type float double)
			set(KRN_NAME "${_cl_kernel}-${_type}")
			set(KRN_VNAME "${_cl_kernel}_${_type}")
			set(KRN_FNAME "src/cl/${KRN_NAME}.cl")
			set(KRN_FULL_FNAME "${CMAKE_CURRENT_SOURCE_DIR}/src/cl/${KRN_NAME}.cl")
			set(KRN_HEADER_FNAME "${CMAKE_CURRENT_BINARY_DIR}/profit/cl/${KRN_NAME}.h")
			add_custom_command(OUTPUT ${KRN_HEADER_FNAME}
				COMMAND ${CMAKE_COMMAND}
				ARGS -DKRN_NAME=${KRN_NAME} -DKRN_VNAME=${KRN_VNAME} -DKRN_FNAME=${KRN_FNAME}
			         -DKRN_FULL_FNAME=${KRN_FULL_FNAME} -DKRN_HEADER_FNAME=${KRN_HEADER_FNAME}
			         -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/GenerateOpenCLHeader.cmake
			    DEPENDS ${KRN_FULL_FNAME}
			)
			list(APPEND OPENCL_KERNEL_HEADERS ${KRN_HEADER_FNAME})
		endforeach()
	endforeach()
endif()

#
# The shared library sources
#
set(PROFIT_SRC
   src/brokenexponential.cpp
   src/convolve.cpp
   src/coresersic.cpp
   src/crc.cpp
   src/exceptions.cpp
   src/ferrer.cpp
   src/fft.cpp
   src/image.cpp
   src/library.cpp
   src/king.cpp
   src/model.cpp
   src/moffat.cpp
   src/opencl.cpp
   src/profile.cpp
   src/psf.cpp
   src/radial.cpp
   src/sersic.cpp
   src/sky.cpp
   src/utils.cpp
)
set(LIB_TYPE SHARED)
if (MSVC)
	set(LIB_TYPE STATIC)
endif()
add_library(profit ${LIB_TYPE} ${PROFIT_SRC} ${OPENCL_KERNEL_HEADERS})
target_link_libraries(profit ${profit_LIBS})
set(PROFIT_TARGETS profit)

#
# The executable
#
add_executable(profit-cli src/fits_utils.cpp src/profit-cli.cpp)
target_link_libraries(profit-cli profit)
set(PROFIT_TARGETS ${PROFIT_TARGETS} profit-cli)

#
# Installing lib + binary + headers
#
install(TARGETS ${PROFIT_TARGETS}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install(FILES
       ${CMAKE_CURRENT_BINARY_DIR}/profit/config.h
        include/profit/common.h
        include/profit/convolve.h
        include/profit/exceptions.h
        include/profit/fft.h
        include/profit/image.h
        include/profit/library.h
        include/profit/model.h
        include/profit/opencl.h
        include/profit/profile.h
        include/profit/profit.h
        include/profit/utils.h
        DESTINATION include/profit)

#
# Unit tests
#
if( LIBPROFIT_TEST )
	find_package(CxxTest)

	# In cmake < 3.7 the cxxtestgen executable is always executed through an
	# interpreter, and doesn't honour the shebang line that could be present.
	# We try to rectify this situation here with a solution similar to what
	# is found in cmake >= 3.7
	if (${CMAKE_MAJOR_VERSION} EQUAL 3 AND ${CMAKE_MINOR_VERSION} LESS 7 AND
	    NOT "${CXXTEST_PYTHON_TESTGEN_EXECUTABLE}" STREQUAL "")
		execute_process(COMMAND ${CXXTEST_PYTHON_TESTGEN_EXECUTABLE} --version
		                OUTPUT_VARIABLE _CXXTEST_OUT ERROR_VARIABLE _CXXTEST_OUT
		                RESULT_VARIABLE _CXXTEST_RESULT)
		if (_CXXTEST_RESULT EQUAL 0)
			set(CXXTEST_TESTGEN_INTERPRETER "")
		endif()
	endif()


	if( CXXTEST_FOUND )
		include(CTest)
		enable_testing()
		add_subdirectory(tests)

		add_test("cli-version" profit-cli -V)
		add_test("cli-help" profit-cli -h)
		add_test("cli-null-profile" profit-cli -p null)
		add_test("cli-many-profiles" profit-cli -p 100,null -p 5,sky)
		add_test("cli-unknown-parameter" profit-cli -p null:unknown=1)
		add_test("cli-text-output" profit-cli -p null -t)
		add_test("cli-fits-output" profit-cli -p null -f image)
		add_test("cli-convolution-with-cmdline" profit-cli -p null:convolve=1 -P 2:2:1,1,2,2)
		add_test("cli-complex-usage" profit-cli -p null -w 100 -H 50 -S 2 -F -m 1 -i 2 -u -T brute -e 1 -r)
		# Needs 'cli-fits-output' to be run first
		add_test("cli-convolution-with-file" profit-cli -p null:convolve=1 -P image.fits)
		if (PROFIT_OPENCL)
			add_test("cli-opencl-list" profit-cli -c)
		endif()
	endif()
endif()
