/**
 * OpenCL utility methods for libprofit
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
#include <fstream>
#include <streambuf>
#include <sstream>
#include <string>
#include <vector>
#include <sys/time.h>

#include "profit/exceptions.h"
#include "profit/opencl.h"

using namespace std;

namespace profit {

map<pair<int, string>, map<int, string>> get_opencl_info() {

	vector<cl::Platform> all_platforms;
	if( cl::Platform::get(&all_platforms) != CL_SUCCESS ) {
		throw opencl_error("Error while getting OpenCL platforms");
	}

	map<pair<int, string>, map<int, string>> pinfo;
	unsigned int pidx = 0;
	for(auto platform: all_platforms) {
		vector<cl::Device> devices;
		platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);

		map<int, string> dinfo;
		unsigned int didx = 0;
		for(auto device: devices) {
			dinfo[didx] = device.getInfo<CL_DEVICE_NAME>();
		}

		string name = platform.getInfo<CL_PLATFORM_NAME>();
		pinfo[make_pair(pidx++, name)] = dinfo;
	}

	return pinfo;
}

static
shared_ptr<OpenCL_env> _get_opencl_environment(unsigned int platform_idx, unsigned int device_idx, bool use_double) {

	vector<cl::Platform> all_platforms;
	if( cl::Platform::get(&all_platforms) != CL_SUCCESS ) {
		throw opencl_error("Error while getting OpenCL platforms");
	}
	if( all_platforms.size() == 0 ){
		throw opencl_error("No platforms found. Check OpenCL installation");
	}

	if( platform_idx >= all_platforms.size() ) {
		ostringstream ss;
		ss << "OpenCL platform index " << platform_idx << " must be < " << all_platforms.size();
		throw invalid_parameter(ss.str());
	}

	cl::Platform platform = all_platforms[platform_idx];

	//get default device of the default platform
	vector<cl::Device> all_devices;
	platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
	if( all_devices.size() == 0 ){
		throw opencl_error("No devices found. Check OpenCL installation");
	}
	if( device_idx >= all_devices.size() ) {
		ostringstream ss;
		ss << "OpenCL device index " << device_idx << " must be < " << all_devices.size();
		throw invalid_parameter(ss.str());
	}

	cl::Device device = all_devices[device_idx];

	if( use_double ) {
		auto config = device.getInfo<CL_DEVICE_DOUBLE_FP_CONFIG>();
		if( config == 0 ) {
			throw opencl_error("Double precision requested but not supported by device");
		}
	}

	// kernel calculates sersic profile for each stuff
	const char *sersic_float =
#include "profit/cl/sersic-float.cl"
	;
	const char *sersic_double =
#include "profit/cl/sersic-double.cl"
	;

	cl::Program::Sources sources;
	sources.push_back(sersic_float);
	if( use_double ) {
		sources.push_back(sersic_double);
	}

	cl::Context context(device);
	cl::Program program(context, sources);
	try {
		program.build({device});
	} catch (const cl::Error &e) {
		throw opencl_error("Error building program: " + program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device));
	}

	cl::CommandQueue queue(context, device);

	return make_shared<OpenCL_env>(OpenCL_env{context, device, queue, program, use_double});
}

shared_ptr<OpenCL_env> get_opencl_environment(unsigned int platform_idx, unsigned int device_idx, bool use_double) {

	// Wrap cl::Error exceptions
	try {
		return _get_opencl_environment(platform_idx, device_idx, use_double);
	} catch(const cl::Error &e) {
		ostringstream os;
		os << "OpenCL error: " << e.what() << ". OpenCL error code: " << e.err();
		throw opencl_error(os.str());
	}
}

}

#endif /* PROFIT_OPENCL */