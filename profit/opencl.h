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

#ifndef PROFIT_OPENCL_H
#define PROFIT_OPENCL_H

#ifdef PROFIT_OPENCL

#include <map>
#include <string>

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

OpenCL_env *get_opencl_environment(unsigned int platform_idx, unsigned int device_idx, bool use_double);

void free_opencl_environment(OpenCL_env *env);

} /* namespace profit */

#endif /* PROFIT_OPENCL */

#endif /* PROFIT_MODEL_H */