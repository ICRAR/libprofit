/**
 * Implementation of library-related routines for libprofit
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

#include <iostream>

#include "profit/config.h"
#include "profit/library.h"

#ifdef PROFIT_FFTW
#include <fftw3.h>
#endif // PROFIT_FFTW

namespace profit {

static std::string _version = std::to_string(PROFIT_VERSION_MAJOR) + "." +
                              std::to_string(PROFIT_VERSION_MINOR) + "." +
                              std::to_string(PROFIT_VERSION_PATCH);

std::string version()
{
	return _version;
}

unsigned short version_major()
{
	return PROFIT_VERSION_MAJOR;
}

unsigned short version_minor()
{
	return PROFIT_VERSION_MINOR;
}

unsigned short version_patch()
{
	return PROFIT_VERSION_PATCH;
}


bool init()
{

	// Initialize FFTW library, including its OpenMP support
#ifdef PROFIT_FFTW
	fftw_import_system_wisdom();
#ifdef PROFIT_FFTW_OPENMP
	int res = fftw_init_threads();
	if (!res) {
		std::cerr << "Error while initializing FFTW threads support, errno = " << res << std::endl;
		return false;
	}
#endif // PROFIT_FFTW_OPENMP
#endif // PROFIT_FFTW

	return true;
}

void finish()
{
#ifdef PROFIT_FFTW
#ifdef PROFIT_FFTW_OPENMP
	fftw_cleanup_threads();
#endif /* PROFIT_FFTW_OPENMP */
	fftw_cleanup();
#endif // PROFIT_FFTW
}

bool has_openmp()
{
#ifdef PROFIT_OPENMP
	return true;
#else
	return false;
#endif // PROFIT_OPENMP
}

bool has_fftw()
{
#ifdef PROFIT_FFTW
	return true;
#else
	return false;
#endif // PROFIT_FFTW
}

bool has_opencl()
{
#ifdef PROFIT_OPENCL
	return true;
#else
	return false;
#endif // PROFIT_OPENCL
}

unsigned short opencl_version_major()
{
#ifdef PROFIT_OPENCL
	return PROFIT_OPENCL_MAJOR;
#else
	return 0;
#endif // PROFIT_OPENCL
}

unsigned short opencl_version_minor()
{
#ifdef PROFIT_OPENCL
	return PROFIT_OPENCL_MINOR;
#else
	return 0;
#endif // PROFIT_OPENCL
}

}