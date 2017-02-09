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
#  define CL_HPP_TARGET_OPENCL_VERSION  120
#  define CL_HPP_MINIMUM_OPENCL_VERSION 110
# endif
# include <CL/cl2.hpp>
#else
# error "We currently support only the OpenCL 2.0 C++ API"
#endif /* PROFIT_OPENCL_MAJOR_VERSION */

namespace profit
{

/**
 * A datatype for storing an OpenCL version.
 * It should have the form major*100 + minor*10 (e.g., 120 for OpenCL 1.2)
 */
typedef unsigned int cl_ver_t;

/**
 * An OpenCL environment
 *
 * This structure holds all the required information to make libprofit work
 * against a given device in a particular platform.
 */
typedef struct _OpenCL_env {

	/** The device to be used throughout OpenCL operations */
	cl::Device device;

	/** The OpenCL supported by the platform this device belongs to */
	cl_ver_t version;

	/** The OpenCL context used throughout the OpenCL operations */
	cl::Context context;

	/** The queue set up against this device to be used by libprofit */
	cl::CommandQueue queue;

	/**
	 * The set of kernels and routines compiled against this device and
	 * required by libprofit
	 */
	cl::Program program;

	/**
	 * Whether double floating-point precision has been requested on this device
	 * or not.
	 */
	bool use_double;

} OpenCL_env;

/**
 * A structure holding information about a specific OpenCL device
 */
typedef struct _OpenCL_dev_info {

	/** The name of the device */
	std::string name;

	/** Whether or not this device supports double floating-point precision */
	bool double_support;

} OpenCL_dev_info;

/**
 * An structure holding information about a specific OpenCL platform.
 */
typedef struct _OpenCL_plat_info {

	/** The name of the platform */
	std::string name;

	/** The supported OpenCL version */
	cl_ver_t supported_opencl_version;

	/** A map containing information about all devices on this platform */

	std::map<int, OpenCL_dev_info> dev_info;
} OpenCL_plat_info;

/**
 * Queries the system about the OpenCL supported platforms and devices and returns
 * the information the caller.
 *
 * @return A map keyed by index, containing the information of each of the
 *         OpenCL platforms found on this system.
 */
std::map<int, OpenCL_plat_info> get_opencl_info();

/**
 * Prepares an OpenCL working space for using with libprofit.
 *
 * This method will get the requested device on the requested platform, compile
 * the libprofit OpenCL kernel sources to be used against it, and set up a queue
 * on the device.
 *
 * @param platform_idx The index of the platform to use
 * @param device_idx The index of device to use in the platform
 * @param use_double Whether double floating-point support should be used in
 *        the device or not.
 * @return A pointer to a OpenCL_env structure, which contains the whole set of
 *         elements required to work with the requested device.
 */
std::shared_ptr<OpenCL_env> get_opencl_environment(unsigned int platform_idx, unsigned int device_idx, bool use_double);

} /* namespace profit */

#endif /* PROFIT_MODEL_H */

#endif /* PROFIT_OPENCL */
