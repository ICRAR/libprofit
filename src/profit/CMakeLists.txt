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


#
# Create the final configuration file that records
# the compile-time options used for this libprofit build
#
configure_file(config.h.in config.h @ONLY)

set(PUBLIC_HEADERS
    ${CMAKE_CURRENT_BINARY_DIR}/config.h
    common.h
    convolve.h
    exceptions.h
    fft.h
    image.h
    library.h
    model.h
    opencl.h
    profile.h
    profit.h
    utils.h
)

set(PRIVATE_HEADERS
    brokenexponential.h
    convolver_impl.h
    coresersic.h
    crc.h
    dot_product.h
    ferrer.h
    fft_impl.h
    fits_utils.h
    king.h
    moffat.h
    null.h
    omp_utils.h
    opencl_impl.h
    psf.h
    radial.h
    sersic.h
    sky.h
)

set(SRC
    brokenexponential.cpp
    convolve.cpp
    coresersic.cpp
    crc.cpp
    exceptions.cpp
    ferrer.cpp
    fft.cpp
    image.cpp
    library.cpp
    king.cpp
    model.cpp
    moffat.cpp
    opencl.cpp
    profile.cpp
    psf.cpp
    radial.cpp
    sersic.cpp
    sky.cpp
    utils.cpp
)

#
# Generate header files with the OpenCL kernel sources out of each individual
# .cl file
#
if (PROFIT_OPENCL)
	foreach(_cl_kernel brokenexponential common convolve coresersic ferrer king moffat sersic)
		foreach(_type float double)
			set(KRN_NAME "${_cl_kernel}-${_type}")
			set(KRN_VNAME "${_cl_kernel}_${_type}")
			set(KRN_FNAME "cl/${KRN_NAME}.cl")
			set(KRN_FULL_FNAME "${CMAKE_CURRENT_SOURCE_DIR}/cl/${KRN_NAME}.cl")
			set(KRN_HEADER_FNAME "${CMAKE_CURRENT_BINARY_DIR}/cl/${KRN_NAME}.h")
			add_custom_command(OUTPUT ${KRN_HEADER_FNAME}
				COMMAND ${CMAKE_COMMAND}
				ARGS -DKRN_NAME=${KRN_NAME} -DKRN_VNAME=${KRN_VNAME} -DKRN_FNAME=${KRN_FNAME}
				      -DKRN_FULL_FNAME=${KRN_FULL_FNAME} -DKRN_HEADER_FNAME=${KRN_HEADER_FNAME}
				      -P ${PROJECT_SOURCE_DIR}/cmake/GenerateOpenCLHeader.cmake
				DEPENDS ${KRN_FULL_FNAME}
			)
		 list(APPEND PRIVATE_HEADERS ${KRN_HEADER_FNAME})
		endforeach()
	endforeach()
endif()

set(LIB_TYPE SHARED)
if (MSVC)
	set(LIB_TYPE STATIC)
endif()
add_library(profit ${LIB_TYPE} ${SRC} ${PUBLIC_HEADERS} ${PRIVATE_HEADERS})
target_link_libraries(profit PRIVATE ${PROFIT_LIBS})
target_include_directories(profit PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_include_directories(profit PUBLIC ${PROJECT_BINARY_DIR}/src)
target_include_directories(profit PRIVATE ${PROFIT_INCLUDE_DIRS})


#
# Installing lib and headers
#
install(TARGETS profit
    EXPORT profit-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(FILES
    ${PUBLIC_HEADERS}
    DESTINATION include/profit
)


if (LIBPROFIT_TEST)
	add_subdirectory(tests)
endif()
