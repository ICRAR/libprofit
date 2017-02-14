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
#include <map>
#include <sstream>
#include <tuple>
#include <vector>

#include "profit/common.h"
#include "profit/exceptions.h"
#include "profit/opencl.h"
#include "profit/model.h"
#include "profit/radial.h"
#include "profit/utils.h"

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

	stats = shared_ptr<RadialProfileStats>(new RadialProfileStats());
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
	auto env = model.opencl_env;
	if( !env || !supports_opencl() ) {
		evaluate_cpu(image);
		return;
	}

	try {
		if( env->use_double ) {
			evaluate_opencl<double>(image);
		}
		else {
			evaluate_opencl<float>(image);
		}
	} catch (const cl::Error &e) {
		ostringstream os;
		os << "OpenCL error: " << e.what() << ". OpenCL error code: " << e.err();
		throw opencl_error(os.str());
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

/*
 * Small trait that describes specific floating types
 */
template <typename T>
struct float_traits {
	const static bool is_float = false;
	const static bool is_double = false;
	constexpr const static char * name = "unknown";
};
template <>
struct float_traits<float> {
	const static bool is_float = true;
	const static bool is_double = false;
	constexpr const static char * name = "float";
};
template <>
struct float_traits<double> {
	const static bool is_float = false;
	const static bool is_double = true;
	constexpr const static char * name = "double";
};

/*
 * Simple structure holding a 2D point
 */
template <typename FT>
struct point_t {
	FT x;
	FT y;
};

/*
 * A structure to hold the information needed to perform subsampling
 * on a specific point
 */
template <typename FT>
struct ss_info_t {
	point_t<FT> point;
	FT xbin;
	FT ybin;
	unsigned int resolution;
	unsigned int max_recursion;
};

/*
 * A smaller version of the previous structure
 * This is what gets passed down to the kernel
 * (and thus we do less memory alloc/copy on the device)
 */
template <typename FT>
struct ss_kinfo_t {
	point_t<FT> point;
	FT xbin;
	FT ybin;
};

template <typename FT>
static inline
unsigned int new_subsampling_points(const vector<ss_info_t<FT>> &last_ss_info, vector<ss_info_t<FT>> &ss_info, unsigned int recur_level) {

	ss_info.clear();

	unsigned int subsampled_pixels = 0;
	for(const auto &info: last_ss_info) {

		const unsigned int res = info.resolution;
		const unsigned int maxr = info.max_recursion;
		FT x = info.point.x;
		FT y = info.point.y;

		if( x == -1 || recur_level >= info.max_recursion) {
			continue;
		}

		FT x0 = x - info.xbin/2;
		FT y0 = y - info.ybin/2;
		FT ss_res_x = info.xbin / res;
		FT ss_res_y = info.ybin / res;

		// we can't cope with more subsampling, sorry
		if( ss_res_x == 0 || ss_res_y == 0 ) {
			continue;
		}

		subsampled_pixels++;
		for(unsigned int j=0; j!=res; j++) {
			FT y_diff = (j*res + 0.5f) * ss_res_y;
			for(unsigned int i=0; i!=res; i++) {
				FT x_diff = (i + 0.5f)*ss_res_x;
				ss_info.push_back({
					{x0 + x_diff, y0 + y_diff},
					ss_res_x, ss_res_y,
					res, maxr
				});
			}
		}
	}

	return subsampled_pixels;
}

template <typename FT>
void RadialProfile::evaluate_opencl(vector<double> &image) {

#define AS_FT(x) static_cast<FT>(x)

	typedef point_t<FT> point_t;
	typedef ss_info_t<FT> ss_info_t;
	typedef ss_kinfo_t<FT> ss_kinfo_t;

	unsigned int imsize = model.width * model.height;

	auto env = model.opencl_env;

	/* Prepare the initial evaluation kernel */
	unsigned int arg = 0;
	auto kname = name + "_" + float_traits<FT>::name;
	cl::Buffer image_buffer(env->context, CL_MEM_WRITE_ONLY, sizeof(FT)*imsize);
	cl::Buffer subsampling_points_buffer(env->context, CL_MEM_WRITE_ONLY, sizeof(point_t)*imsize);
	cl::Kernel kernel = cl::Kernel(env->program, kname.c_str());
	kernel.setArg(arg++, image_buffer);
	kernel.setArg(arg++, subsampling_points_buffer);
	kernel.setArg(arg++, model.width);
	kernel.setArg(arg++, model.height);
	kernel.setArg(arg++, (int)rough);
	kernel.setArg(arg++, AS_FT(model.scale_x));
	kernel.setArg(arg++, AS_FT(model.scale_y));
	add_common_kernel_parameters<FT>(arg, kernel);

	cl::Event fill_im_evt, fill_ss_points_evt, kernel_evt, read_evt, read_ss_points_evt;

	// OpenCL 1.2 allows to do this; otherwise the work has to be done in the kernel
	// (which we do)
	cl::vector<cl::Event> k_wait_evts;
	if( env->version >= 120 ) {
		env->queue.enqueueFillBuffer<FT>(image_buffer, 0, 0, sizeof(FT)*imsize, NULL, &fill_im_evt);
		env->queue.enqueueFillBuffer<point_t>(subsampling_points_buffer, {-1, -1}, 0, sizeof(point_t)*imsize, NULL, &fill_ss_points_evt);
		k_wait_evts.push_back(fill_im_evt);
		k_wait_evts.push_back(fill_ss_points_evt);
	}

	// Enqueue the kernel, and read back the resulting image + set of points to subsample
	env->queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(imsize), cl::NullRange, &k_wait_evts, &kernel_evt);

	cl::vector<cl::Event> read_waiting_evts{kernel_evt};
	if( float_traits<FT>::is_double ) {
		env->queue.enqueueReadBuffer(image_buffer, CL_FALSE, 0, sizeof(double)*imsize, image.data(), &read_waiting_evts, &read_evt);
		read_evt.wait();
	}
	else {
		vector<FT> image_from_kernel(image.size());
		env->queue.enqueueReadBuffer(image_buffer, CL_FALSE, 0, sizeof(FT)*imsize, image_from_kernel.data(), &read_waiting_evts, &read_evt);
		read_evt.wait();
		copy(image_from_kernel.begin(), image_from_kernel.end(), image.begin());
	}

	// we're done here
	if( rough ) {
		return;
	}

	vector<point_t> ss_points(image.size());
	env->queue.enqueueReadBuffer(subsampling_points_buffer, CL_TRUE, 0, sizeof(point_t)*imsize, ss_points.data(), &read_waiting_evts, &read_ss_points_evt);
	read_ss_points_evt.wait();

	// enrich the points to subsample with their subsampling information
	vector<ss_info_t> last_ss_info;
	last_ss_info.reserve(image.size());
	FT half_xbin = AS_FT(model.scale_x)/2;
	FT half_ybin = AS_FT(model.scale_y)/2;
	unsigned int top_recursions = 0;
	for(auto const &point: ss_points) {
		if( point.x == -1 ) {
			continue;
		}
		unsigned int resolution, max_recursions;
		subsampling_params(point.x, point.y, resolution, max_recursions);
		top_recursions = max(top_recursions, max_recursions);
		last_ss_info.push_back({point, AS_FT(model.scale_x), AS_FT(model.scale_y), resolution, max_recursions});
	}

	auto ss_kname = name + "_subsample_" + float_traits<FT>::name;
	cl::Kernel subsample_kernel = cl::Kernel(env->program, ss_kname.c_str());

	typedef struct _im_result {
		point_t point;
		FT value;
	} im_result_t;
	vector<im_result_t> subimages_results;

	// Preparing for the recursive subsampling
	vector<ss_info_t> ss_info;
	unsigned int recur_level = 0;
	unsigned int total_subsamples = 0;

	while( recur_level < top_recursions ) {

		unsigned int subsampled_pixels = new_subsampling_points<FT>(last_ss_info, ss_info, recur_level);

		auto subsamples = ss_info.size();
		if( !subsamples ) {
			break;
		}

		total_subsamples += subsamples;

#ifdef PROFIT_DEBUG
		/* record how many sub-integrations we've done */
		n_integrations[recur_level] = subsampled_pixels;
#endif

		/* Keeping things in size */
		auto last_im_idx = subimages_results.size();
		subimages_results.reserve(last_im_idx + subsamples);
		last_ss_info.resize(subsamples);

		try {

			cl::Buffer subimage_buf(env->context, CL_MEM_WRITE_ONLY, sizeof(FT)*subsamples);
			cl::Buffer ss_kinfo_buf(env->context, CL_MEM_READ_WRITE, sizeof(ss_kinfo_t)*subsamples);

			arg = 0;
			subsample_kernel.setArg(arg++, subimage_buf);
			subsample_kernel.setArg(arg++, ss_kinfo_buf);
			subsample_kernel.setArg(arg++, AS_FT(acc));
			add_common_kernel_parameters<FT>(arg, subsample_kernel);

			cl::Event kernel_evt, w_ss_kinfo_evt, r_kimage_evt, r_ss_kinfo_evt, fill_evt;
			vector<FT> kimage(subsamples);

			// The information we pass down to the kernels is a subset of the original
			vector<ss_kinfo_t> ss_kinfo(subsamples);
			transform(ss_info.begin(), ss_info.end(), ss_kinfo.begin(), [](const ss_info_t &info) {
				return ss_kinfo_t{info.point, info.xbin, info.ybin};
			});

			if( env->version >= 120 ) {
				env->queue.enqueueFillBuffer<FT>(subimage_buf, 0, 0, sizeof(FT)*subsamples, NULL, &fill_evt);
			}
			env->queue.enqueueWriteBuffer(ss_kinfo_buf, CL_FALSE, 0, sizeof(ss_kinfo_t)*subsamples, ss_kinfo.data(), NULL, &w_ss_kinfo_evt);
			cl::vector<cl::Event> kernel_waiting_evts{w_ss_kinfo_evt};
			if( env->version >= 120 ) {
				kernel_waiting_evts.push_back(fill_evt);
			}
			env->queue.enqueueNDRangeKernel(subsample_kernel, cl::NullRange, cl::NDRange(subsamples), cl::NullRange, &kernel_waiting_evts, &kernel_evt);
			cl::vector<cl::Event> read_waiting_evts{kernel_evt};
			env->queue.enqueueReadBuffer(subimage_buf, CL_FALSE, 0, sizeof(FT)*subsamples, kimage.data(), &read_waiting_evts, &r_kimage_evt);
			env->queue.enqueueReadBuffer(ss_kinfo_buf, CL_FALSE, 0, sizeof(ss_kinfo_t)*subsamples, ss_kinfo.data(), &read_waiting_evts, &r_ss_kinfo_evt);
			env->queue.finish();

			// Feed back the kinfo to the main subsampling info vectors
			auto ss_info_it = ss_info.begin();
			auto last_ss_info_it  = last_ss_info.begin();
			auto kim_it = kimage.begin();
			subimages_results.reserve(last_im_idx + subsamples);

			for(auto &kinfo: ss_kinfo) {

				// Copy the point information from the kernel
				last_ss_info_it->point = kinfo.point;

				// ... and the rest of the subsampling info for next round
				last_ss_info_it->xbin = ss_info_it->xbin;
				last_ss_info_it->ybin = ss_info_it->ybin;
				last_ss_info_it->resolution = ss_info_it->resolution;
				last_ss_info_it->max_recursion = ss_info_it->max_recursion;

				FT val = *kim_it;
				for(int i=0; i<=recur_level; i++) {
					val /= (ss_info_it->resolution * ss_info_it->resolution);
				}

				// This is a final result that requires no more subsampling
//				if( kinfo.point.x == -1 ) {
//					FT x = ss_info_it->point.x / model.scale_x;
//					FT y = ss_info_it->point.y / model.scale_y;
//					unsigned int idx = static_cast<unsigned int>(floor(x)) + static_cast<unsigned int>(floor(y)) * model.width;
//					image[idx] += val;
//				}
				subimages_results.push_back(im_result_t{ss_info_it->point, val});

				last_ss_info_it++;
				ss_info_it++;
				kim_it++;
			}

		} catch(const cl::Error &e) {
			// running out of memory, cannot go any further
			if( e.err() == CL_INVALID_BUFFER_SIZE ) {
				break;
			}
			throw e;
		}

		recur_level++;
	}

	for_each(subimages_results.begin(), subimages_results.end(), [&image, this](const im_result_t &res) {
		FT x = res.point.x / model.scale_x;
		FT y = res.point.y / model.scale_y;
		unsigned int idx = static_cast<unsigned int>(floor(x)) + static_cast<unsigned int>(floor(y)) * model.width;
		image[idx] += res.value;
	});

	// the image needs to be multiplied by the pixel scale
	double scale = this->get_pixel_scale();
	transform(image.begin(), image.end(), image.begin(), [scale](double pixel) {
		return pixel * scale;
	});

}

template <typename FT>
void RadialProfile::add_common_kernel_parameters(unsigned int arg, cl::Kernel &kernel) const {
	kernel.setArg(arg++, static_cast<FT>(xcen));
	kernel.setArg(arg++, static_cast<FT>(ycen));
	kernel.setArg(arg++, static_cast<FT>(_cos_ang));
	kernel.setArg(arg++, static_cast<FT>(_sin_ang));
	kernel.setArg(arg++, static_cast<FT>(axrat));
	kernel.setArg(arg++, static_cast<FT>(rscale));
	kernel.setArg(arg++, static_cast<FT>(rscale_switch));
	kernel.setArg(arg++, static_cast<FT>(rscale_max));
	kernel.setArg(arg++, static_cast<FT>(box));
	if( float_traits<FT>::is_float ) {
		add_kernel_parameters_float(arg, kernel);
	}
	else {
		add_kernel_parameters_double(arg, kernel);
	}
}

bool RadialProfile::supports_opencl() const {
	return false;
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
void RadialProfile::add_kernel_parameters_float(unsigned int index, cl::Kernel &kernel) const {
	return;
}

void RadialProfile::add_kernel_parameters_double(unsigned int index, cl::Kernel &kernel) const {
	return;
}
#endif /* PROFIT_OPENCL */

} /* namespace profit */