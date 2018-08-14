/**
 * Base Profile class implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
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

#include <sstream>
#include <string>

#include "profit/common.h"
#include "profit/exceptions.h"
#include "profit/profile.h"


namespace profit {

ProfileStats::ProfileStats() :
	total(0)
{
	// no-op
}

ProfileStats::~ProfileStats()
{
	// no-op
}

RadialProfileStats::RadialProfileStats() :
	ProfileStats()
#ifdef PROFIT_OPENCL
	,cl_times(),
	subsampling{0, 0, 0, OpenCL_times(), 0, 0},
	final_image(0)
#endif /* PROFIT_OPENCL */
{
	// no-op
}

Profile::Profile(const Model &model, const std::string &name) :
	model(model),
	name(name),
	stats(),
	convolve(false)
{
	register_parameter("convolve", convolve);
}

Profile::~Profile()
{
	// no-op
}

bool Profile::do_convolve() const {
	return convolve;
}

const std::string& Profile::get_name() const {
	return name;
}

std::shared_ptr<ProfileStats> Profile::get_stats() const {
	return stats;
}

void Profile::register_parameter(const char *name, bool &parameter)
{
	bool_parameters.insert({name, parameter});
}

void Profile::register_parameter(const char *name, unsigned int &parameter)
{
	uint_parameters.insert({name, parameter});
}

void Profile::register_parameter(const char *name, double &parameter)
{
	double_parameters.insert({name, parameter});
}

template <typename T>
void set_parameter(
	Profile::parameter_holder<T> &parameters,
	const std::string &name,
	const std::string &profile_name,
	T val)
{
	constexpr auto tname = type_info<T>::name;
	if (parameters.find(name) == parameters.end()) {
		std::ostringstream os;
		os << "Unknown " << tname << " parameter in profile " << profile_name << ": " << name;
		throw invalid_parameter(os.str());
	}
	parameters.at(name).get() = val;
}

void Profile::parameter(const std::string &name, bool val) {
	set_parameter(bool_parameters, name, get_name(), val);
}

void Profile::parameter(const std::string &name, double val) {
	set_parameter(double_parameters, name, get_name(), val);
}

void Profile::parameter(const std::string &name, unsigned int val) {
	set_parameter(uint_parameters, name, get_name(), val);
}

} /* namespace profit */