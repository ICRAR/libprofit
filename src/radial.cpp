/**
 * Radial profile base implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Aaron Robotham, Dan Taranu, Rodrigo Tobar
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

#include <algorithm>
#include <cmath>
#include <vector>
#include <tuple>

#include "profit/common.h"
#include "profit/exceptions.h"
#include "profit/model.h"
#include "profit/radial.h"
#include "profit/utils.h"

#ifdef PROFIT_OPENCL
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION  120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#include <CL/cl2.hpp>
#include "profit/opencl.h"
#endif /* PROFIT_OPENCL */

using namespace std;

namespace profit
{

inline
void RadialProfile::_image_to_profile_coordinates(double x, double y, double &x_prof, double &y_prof) {
	x -= this->xcen;
	y -= this->ycen;
	x_prof =  x * this->_cos_ang + y * this->_sin_ang;
	y_prof = -x * this->_sin_ang + y * this->_cos_ang;
	y_prof /= this->axrat;
}

double RadialProfile::subsample_pixel(double x0, double x1, double y0, double y1,
                                      unsigned int recur_level, unsigned int max_recursions,
                                      unsigned int resolution) {

	double xbin = (x1-x0) / resolution;
	double ybin = (y1-y0) / resolution;
	double half_xbin = xbin/2.;
	double half_ybin = ybin/2.;
	double total = 0, subval, testval;
	double x , y, x_prof, y_prof;
	unsigned int i, j;

	bool recurse = resolution > 1 && recur_level < max_recursions;

#ifdef PROFIT_DEBUG
	/* record how many sub-integrations we've done */
	if( n_integrations.find(recur_level) != n_integrations.end() ) {
		n_integrations[recur_level] += 1;
	}
	else {
		n_integrations[recur_level] = 1;
	}
#endif

	/* The middle X/Y value is used for each pixel */
	vector<tuple<double, double>> subsample_points;
	x = x0;

	vector<unsigned int> idxs(resolution * resolution);
	if( recurse ) {
		for(i=0; i < resolution; i++) {
			x += half_xbin;
			y = y0;
			for(j=0; j < resolution; j++) {
				y += half_ybin;

				this->_image_to_profile_coordinates(x, y, x_prof, y_prof);
				subval = this->evaluate_at(x_prof, y_prof);

				double delta_y_prof = (-xbin*this->_sin_ang + ybin*this->_cos_ang)/this->axrat;
				testval = this->evaluate_at(abs(x_prof), abs(y_prof) + abs(delta_y_prof));
				if( abs(testval/subval - 1.0) > this->acc ) {
					subsample_points.push_back(make_tuple(x, y));
				}
				else {
					total += subval;
				}
				y += half_ybin;
			}

			x += half_xbin;
		}
	}
	else {
		for(i=0; i < resolution; i++) {
			x += half_xbin;
			y = y0;
			for(j=0; j < resolution; j++) {
				y += half_ybin;
				this->_image_to_profile_coordinates(x, y, x_prof, y_prof);
				total += this->evaluate_at(x_prof, y_prof);
				y += half_ybin;
			}
			x += half_xbin;
		}
	}

	for(auto &point: subsample_points) {
		double x = get<0>(point);
		double y = get<1>(point);
		total += this->subsample_pixel(x - half_xbin, x + half_xbin,
		                               y - half_ybin, y + half_ybin,
		                               recur_level + 1, max_recursions,
		                               resolution);
	}

	/* Average and return */
	return total / (resolution * resolution);
}

void RadialProfile::initial_calculations() {

	/*
	 * get_rscale() is implemented by subclasses. It provides the translation
	 * from profile-specific parameters into the common rscale concept used in
	 * this common class.
	 */
	this->rscale = this->get_rscale();

	/*
	 * Calculate the total luminosity used by this profile, used
	 * later to calculate the exact contribution of each pixel.
	 */
	double box = this->box + 2;
	double r_box = M_PI * box / (2*beta(1/box, 1/box));
	double lumtot = this->get_lumtot(r_box);
	this->_ie = pow(10, -0.4*(this->mag - this->model.magzero))/lumtot;

	/*
	 * Optionally adjust the user-given rscale_switch and resolution parameters
	 * to more sensible values that will result in faster profile calculations.
	 */
	if( this->adjust ) {

		/*
		 * Automatially adjust the rscale_switch.
		 * Different profiles do it in different ways
		 */
		this->rscale_switch = this->adjust_rscale_switch();

		/*
		 * Calculate a bound, adaptive upscale
		 */
		unsigned int resolution;
		resolution = (unsigned int)ceil(160 / (this->rscale_switch * this->rscale));
		resolution += resolution % 2;
		resolution = max(4, min(16, (int)resolution));
		this->resolution = resolution;

		/*
		 * If the user didn't give a rscale_max we calculate one that covers
		 * %99.99 of the flux
		 */
		if( this->rscale_max == 0 ) {
			this->rscale_max = this->adjust_rscale_max();
		}

		/* Adjust the accuracy we'll use for sub-pixel integration */
		this->acc = this->adjust_acc();

	}

	/*
	 * Get the rotation angle in radians and calculate the coefficients
	 * that will fill the rotation matrix we'll use later to transform
	 * from image coordinates into profile coordinates.
	 *
	 * In galfit the angle started from the Y image axis.
	 */
	double angrad = fmod(this->ang + 90, 360.) * M_PI / 180.;
	this->_cos_ang = cos(angrad);
	this->_sin_ang = sin(angrad);

}

/**
 * The profile validation function
 */
void RadialProfile::validate() {
	if ( axrat <= 0 ) {
		throw invalid_parameter("axrat <= 0, must have axrat > 0");
	}
	if ( axrat > 1 ) {
		throw invalid_parameter("axrat > 1, must have axrat <= 1");
	}
	if ( box <= -2 ) {
		throw invalid_parameter("box <= -2, must have box > -2");
	}
}

/**
 * The scale by which each image pixel value is multiplied
 */
double RadialProfile::get_pixel_scale() {
	double pixel_area = this->model.scale_x * this->model.scale_y;
	return pixel_area * this->_ie;
}

void RadialProfile::subsampling_params(double x, double y,
                                       unsigned int &resolution,
                                       unsigned int &max_recursions) {
	resolution = this->resolution;
	max_recursions = this->max_recursions;
}

/**
 * The main profile evaluation function
 */
void RadialProfile::evaluate(vector<double> &image) {

	/*
	 * Perform all the pre-calculations needed by the radial profiles
	 * (e.g., Ie, cos/sin ang, etc).
	 * We store these profile-global results in the profile object itself
	 * (it contains extra members to store these values) to avoid passing a long
	 * list of values around every method call.
	 */
	this->initial_calculations();

#ifdef PROFIT_DEBUG
	n_integrations.clear();
#endif

#ifndef PROFIT_OPENCL
	evaluate_cpu(image);
#else
	/*
	 * We fallback to the CPU implementation if no OpenCL context has been
	 * given, or if there is no OpenCL kernel implementing the profile
	 */
	OpenCL_env *env = model.opencl_env;
	if( !env ) {
		evaluate_cpu(image);
	}
	else {

		const char *kernel_name;
		if( env->use_double ) {
			kernel_name = get_opencl_kernel_name_double();
		}
		else {
			kernel_name = get_opencl_kernel_name_float();
		}

		if( strlen(kernel_name) == 0 ) {
			evaluate_cpu(image);
		}

		try {
			if( env->use_double ) {
				evaluate_opencl<double>(image, kernel_name);
			}
			else {
				evaluate_opencl<float>(image, kernel_name);
			}
		} catch (const cl::Error &e) {
			throw opencl_error(e.what());
		}
	}
#endif /* PROFIT_OPENCL */

}

void RadialProfile::evaluate_cpu(vector<double> &image) {

	unsigned int i, j;
	double x, y, pixel_val;
	double x_prof, y_prof, r_prof;
	double half_xbin = model.scale_x/2.;
	double half_ybin = model.scale_x/2.;

	double scale = this->get_pixel_scale();

	/* The middle X/Y value is used for each pixel */
	y = 0;
	for(j=0; j < model.height; j++) {
		y += half_ybin;
		x = 0;
		for(i=0; i < model.width; i++) {
			x += half_xbin;

			/* We were instructed to ignore this pixel */
			if( !model.calcmask.empty() && !model.calcmask[i + j*model.width] ) {
				x += half_xbin;
				continue;
			}

			this->_image_to_profile_coordinates(x, y, x_prof, y_prof);

			/*
			 * Check whether we need further refinement.
			 * TODO: the radius calculation doesn't take into account boxing
			 */
			r_prof = sqrt(x_prof*x_prof + y_prof*y_prof);
			if( this->rscale_max > 0 && r_prof/this->rscale > this->rscale_max ) {
				pixel_val = 0.;
			}
			else if( this->rough || r_prof/this->rscale > this->rscale_switch ) {
				pixel_val = this->evaluate_at(x_prof, y_prof);
			}
			else {

				unsigned int resolution;
				unsigned int max_recursions;
				this->subsampling_params(x, y, resolution, max_recursions);

				/* Subsample and integrate */
				pixel_val =  this->subsample_pixel(x - half_xbin, x + half_xbin,
				                                   y - half_ybin, y + half_ybin,
				                                   0, max_recursions, resolution);
			}

			image[i + j*model.width] = scale * pixel_val;
			x += half_xbin;
		}
		y += half_ybin;
	}

}

#ifdef PROFIT_OPENCL

template <typename FT>
void RadialProfile::evaluate_opencl(vector<double> &image, const char *kernel_name) {

	unsigned int i, j;
	double x, y, pixel_val;
	double x_prof, y_prof, r_prof;
	double half_xbin = model.scale_x/2.;
	double half_ybin = model.scale_x/2.;
	unsigned int imsize = model.width * model.height;

	OpenCL_env *env = model.opencl_env;
	double scale = this->get_pixel_scale();

	cl::Buffer buffer_image(env->context, CL_MEM_WRITE_ONLY, sizeof(FT)*imsize);
	cl::Kernel kernel = cl::Kernel(env->program, kernel_name);
	kernel.setArg(0,  buffer_image);
	kernel.setArg(1,  model.width);
	kernel.setArg(2,  model.height);
	kernel.setArg(3,  (FT)xcen);
	kernel.setArg(4,  (FT)ycen);
	kernel.setArg(5,  (FT)_cos_ang);
	kernel.setArg(6,  (FT)_sin_ang);
	kernel.setArg(7,  (FT)axrat);
	kernel.setArg(8,  (FT)rscale);
	kernel.setArg(9,  (FT)rscale_switch);
	kernel.setArg(10, (FT)rscale_max);
	kernel.setArg(11, (int)rough);
	kernel.setArg(12, (FT)box);
	kernel.setArg(13, (FT)scale);
	add_kernel_parameters_float(14, kernel);

	cl::Event kernel_evt;
	env->queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(imsize), cl::NullRange, 0, &kernel_evt);
	read_image_from_kernel<FT>(image, kernel_evt, buffer_image);
}

template <typename FT>
void RadialProfile::read_image_from_kernel(vector<double> &image, cl::Event &kernel_evt, cl::Buffer &buffer_image) const {
	vector<FT> image_from_kernel(image.size());
	cl::vector<cl::Event> read_waiting_evts{kernel_evt};
	model.opencl_env->queue.enqueueReadBuffer(buffer_image, CL_TRUE, 0, sizeof(FT)*image.size(), image_from_kernel.data(), &read_waiting_evts, NULL);
	copy(image_from_kernel.begin(), image_from_kernel.end(), image.begin());
}

template <>
void RadialProfile::read_image_from_kernel<double>(vector<double> &image, cl::Event &kernel_evt, cl::Buffer &buffer_image) const {
	cl::vector<cl::Event> read_waiting_evts{kernel_evt};
	model.opencl_env->queue.enqueueReadBuffer(buffer_image, CL_TRUE, 0, sizeof(double)*image.size(), image.data(), &read_waiting_evts, NULL);
}

#endif /* PROFIT_OPENCL */

/**
 * Constructor with sane defaults
 */
RadialProfile::RadialProfile(const Model &model, const string &name) :
	Profile(model, name),
	xcen(0), ycen(0),
	mag(15), ang(0),
	axrat(1), box(0),
	rough(false), acc(0.1),
	rscale_switch(1), resolution(9),
	max_recursions(2), adjust(true),
	rscale_max(0)
{
	// no-op
}

#ifdef PROFIT_DEBUG
std::map<int,int> RadialProfile::get_integrations() {
	return n_integrations;
}
#endif

bool RadialProfile::parameter_impl(const string &name, bool value) {

	if( Profile::parameter_impl(name, value) ) {
		return true;
	}

	if( name == "rough" )              { rough = value; }
	else if( name == "adjust" )        { adjust = value; }
	else {
		return false;
	}

	return true;
}

bool RadialProfile::parameter_impl(const string &name, double value) {

	if( Profile::parameter_impl(name, value) ) {
		return true;
	}

	if( name == "xcen" )               { xcen = value; }
	else if( name == "ycen" )          { ycen = value; }
	else if( name == "mag" )           { mag = value; }
	else if( name == "ang" )           { ang = value; }
	else if( name == "axrat" )         { axrat = value; }
	else if( name == "box" )           { box = value; }
	else if( name == "acc" )           { acc = value; }
	else if( name == "rscale_switch" ) { rscale_switch = value; }
	else if( name == "rscale_max" )    { rscale_max = value; }
	else {
		return false;
	}

	return true;
}

bool RadialProfile::parameter_impl(const string &name, unsigned int value) {

	if( Profile::parameter_impl(name, value) ) {
		return true;
	}

	if( name == "max_recursions" )  { max_recursions = value; }
	else if( name == "resolution" ) { resolution = value; }
	else {
		return false;
	}

	return true;
}

#ifdef PROFIT_OPENCL
const char * RadialProfile::get_opencl_kernel_name_float() const {
	return "";
}

const char * RadialProfile::get_opencl_kernel_name_double() const {
	return "";
}

void RadialProfile::add_kernel_parameters_float(unsigned int index, cl::Kernel &kernel) const {
	return;
}

void RadialProfile::add_kernel_parameters_double(unsigned int index, cl::Kernel &kernel) const {
	return;
}
#endif /* PROFIT_OPENCL */

} /* namespace profit */