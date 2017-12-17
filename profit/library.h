/**
 * Definitions of library-related routines for libprofit
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2017
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Rodrigo Tobar
 *
 * This file is part of libprofit.
 *
 * libprofit is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libprofit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libprofit.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PROFIT_INIT_FINI_H_
#define PROFIT_INIT_FINI_H_

#include <string>

namespace profit {

/// Returns the version of this libprofit library
/// @return The version of this libprofit library
std::string version();

/// Returns the major version of this libprofit library
/// @return The major version of this libprofit library
unsigned short version_major();

/// Returns the minor version of this libprofit library
/// @return The minor version of this libprofit library
unsigned short version_minor();

/// Returns the patch version of this libprofit library
/// @return The patch version of this libprofit library
unsigned short version_patch();

/// Initializes all static requirements of libprofit.
/// @return If the initialization was correct
bool init();

/// Finalizes all static requirements of libprofit
void finish();

/// Returns whether libprofit was compiled with OpenMP support
/// @return Whether libprofit was compiled with OpenMP support
bool has_openmp();

/// Returns whether libprofit was compiled with FFTW support
/// @return Whether libprofit was compiled with FFTW support
bool has_fftw();

/// Returns whether libprofit was compiled with OpenCL support
/// @return Whether libprofit was compiled with OpenCL support
bool has_opencl();

/// If OpenCL is supported, returns the major portion of the highest OpenCL
/// platform version libprofit can work against. For example, if libprofit was
/// compiled against a platform supporting OpenCL 2.1, this method returns 2.
/// If OpenCL is not supported, the result is undefined.
/// @return The major highest OpenCL platform version that libprofit can work
/// against.
unsigned short opencl_version_major();

/// If OpenCL is supported, returns the minor portion of the highest OpenCL
/// platform version libprofit can work against. For example, if libprofit was
/// compiled against a platform supporting OpenCL 1.2, this method returns 2.
/// If OpenCL is not supported, the result is undefined.
unsigned short opencl_version_minor();

}

#endif // PROFIT_INIT_FINI_H_