/**
 * Header file for OpenCL functionality
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

#ifdef PROFIT_OPENCL

#ifndef PROFIT_OPENCL_H
#define PROFIT_OPENCL_H

#include <map>
#include <memory>
#include <string>

/*
 * OpenCL 1.x uses a different include file, we don't support it yet
 * Also we only give specific inputs to CL/cl2.hpp when building
 * the library, but not when including this file from an external program
 */
#if (PROFIT_OPENCL_MAJOR_VERSION >= 2)
# if defined(PROFIT_BUILD)
#  define CL_HPP_ENABLE_EXCEPTIONS
#  define CL_HPP_TARGET_OPENCL_VERSION  110
#  define CL_HPP_MINIMUM_OPENCL_VERSION 110
# endif
# include <CL/cl2.hpp>
#else
# error "We currently support only the OpenCL 2.0 C++ API"
#endif /* PROFIT_OPENCL_MAJOR_VERSION */

namespace profit
{

typedef struct _OpenCL_env {
	cl::Context context;
	cl::Device device;
	cl::CommandQueue queue;
	cl::Program program;
	bool use_double;
} OpenCL_env;

std::map<std::pair<int, std::string>, std::map<int, std::string>> get_opencl_info();

std::shared_ptr<OpenCL_env> get_opencl_environment(unsigned int platform_idx, unsigned int device_idx, bool use_double);

} /* namespace profit */

#endif /* PROFIT_MODEL_H */

#endif /* PROFIT_OPENCL */
