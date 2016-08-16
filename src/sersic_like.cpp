/**
 * Sersic-like profile implementations
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

#include <cmath>
#include <algorithm>

#include "sersic_like.h"
#include "utils.h"

using namespace std;

namespace profit
{

void SersicLikeProfile::_image_to_profile_coordinates(double x, double y, double &x_prof, double &y_prof) {
	x -= this->xcen;
	y -= this->ycen;
	x_prof =  x * this->_cos_ang + y * this->_sin_ang;
	y_prof = -x * this->_sin_ang + y * this->_cos_ang;
	y_prof /= this->axrat;
}

double SersicLikeProfile::subsample_pixel(double x0, double x1, double y0, double y1,
                                          unsigned int recur_level, unsigned int max_recursions,
                                          unsigned int resolution) {

	double xbin = (x1-x0) / resolution;
	double ybin = (y1-y0) / resolution;
	double half_xbin = xbin/2.;
	double half_ybin = ybin/2.;
	double total = 0, subval, testval;
	double x , y, x_ser, y_ser;
	unsigned int i, j;

	bool recurse = resolution > 1 && recur_level < max_recursions;

#ifdef PROFIT_DEBUG
	/* record how many sub-integrations we've done */
	if( sp->n_integrations.find(recur_level) != sp->n_integrations.end() ) {
		sp->n_integrations[recur_level] += 1;
	}
	else {
		sp->n_integrations[recur_level] = 1;
	}
#endif

	/* The middle X/Y value is used for each pixel */
	x = x0;
	for(i=0; i < resolution; i++) {
		x += half_xbin;
		y = y0;
		for(j=0; j < resolution; j++) {
			y += half_ybin;

			this->_image_to_profile_coordinates(x, y, x_ser, y_ser);
			subval = this->_eval_function(this, x_ser, y_ser, 0, false);

			if( recurse ) {
				double delta_y_ser = (-xbin*this->_sin_ang + ybin*this->_cos_ang)/this->axrat;
				testval = this->_eval_function(this, abs(x_ser), abs(y_ser) + abs(delta_y_ser), 0, false);
				if( abs(testval/subval - 1.0) > this->acc ) {
					subval = this->subsample_pixel(x - half_xbin, x + half_xbin,
					                               y - half_ybin, y + half_ybin,
					                               recur_level + 1, max_recursions,
					                               resolution);
				}
			}

			total += subval;
			y += half_ybin;
		}

		x += half_xbin;
	}

	/* Average and return */
	return total / (resolution * resolution);
}

double SersicLikeProfile::adjust_acc() {
	return 0.1/axrat;
}

void SersicLikeProfile::initial_calculations() {

	/*
	 * get_re() is implemented by subclasses. It provides the translation from
	 * profile-specific parameters into the common "Re" concept used in this
	 * common class.
	 */
	this->_re = this->get_re();

	/*
	 * Calculate the total luminosity used by this profile, used
	 * later to calculate the exact contribution of each pixel.
	 */
	double box = this->box + 2;
	double r_box = M_PI * box / (4*beta(1/box, 1 + 1/box));
	double lumtot = this->get_lumtot(r_box);
	this->_ie = pow(10, -0.4*(this->mag - this->model->magzero))/lumtot;

	/*
	 * Optionally adjust the user-given re_switch (totally) and resolution
	 * (partially) parameters to more sensible values that will result in faster
	 * profile calculations.
	 */
	if( this->adjust ) {

		/*
		 * Automatially adjust the re_switch.
		 * Different profiles do it in different ways
		 */
		this->re_switch = this->adjust_re_switch();

		/*
		 * Calculate a bound, adaptive upscale
		 */
		unsigned int resolution;
		resolution = (unsigned int)ceil(160 / this->re_switch);
		resolution += resolution % 2;
		resolution = max(4, min(16, (int)resolution));
		this->resolution = resolution;

		/*
		 * If the user didn't give a re_max we calculate one that covers
		 * %99.99 of the flux
		 */
		if( this->re_max == 0 ) {
			this->re_max = this->adjust_re_max();
		}

		/* Adjust the accuracy we'll use for sub-pixel integration */
		this->acc = this->adjust_acc();

	}

	/*
	 * Get the rotation angle in radians and calculate the coefficients
	 * that will fill the rotation matrix we'll use later to transform
	 * from image coordinates into sersic coordinates.
	 *
	 * In galfit the angle started from the Y image axis.
	 */
	double angrad = fmod(this->ang + 90, 360.) * M_PI / 180.;
	this->_cos_ang = cos(angrad);
	this->_sin_ang = sin(angrad);

}

/**
 * The sersic validation function
 */
void SersicLikeProfile::validate() {
	// no-op
}

/**
 * The scale by which each image pixel value is multiplied
 */
double SersicLikeProfile::get_pixel_scale() {
	double pixel_area = this->model->scale_x * this->model->scale_y;
	return pixel_area * this->_ie;
}

void SersicLikeProfile::subsampling_params(double x, double y,
                                           unsigned int &resolution,
                                           unsigned int &max_recursions) {
	resolution = this->resolution;
	max_recursions = this->max_recursions;
}

/**
 * The main sersic evaluation function
 */
void SersicLikeProfile::evaluate(double *image) {

	unsigned int i, j;
	double x, y, pixel_val;
	double x_ser, y_ser, r_ser;
	double half_xbin = model->scale_x/2.;
	double half_ybin = model->scale_x/2.;

	this->_eval_function = this->get_evaluation_function();

	/*
	 * All the pre-calculations needed by the sersic-like profiles (Ie, cos/sin ang, etc)
	 * We store these profile-global results in the profile object itself
	 * (it contains extra members to store these values) to avoid passing a long
	 * list of values around every method call.
	 */
	this->initial_calculations();

	double scale = this->get_pixel_scale();

	/* The middle X/Y value is used for each pixel */
	y = 0;
	for(j=0; j < model->height; j++) {
		y += half_ybin;
		x = 0;
		for(i=0; i < model->width; i++) {
			x += half_xbin;

			/* We were instructed to ignore this pixel */
			if( model->calcmask && !model->calcmask[i + j*model->width] ) {
				x += half_xbin;
				continue;
			}

			this->_image_to_profile_coordinates(x, y, x_ser, y_ser);

			/*
			 * Check whether we need further refinement.
			 * TODO: the radius calculation doesn't take into account boxing
			 */
			r_ser = sqrt(x_ser*x_ser + y_ser*y_ser);
			if( this->re_max > 0 && r_ser/this->_re > this->re_max ) {
				pixel_val = 0.;
			}
			else if( this->rough || r_ser/this->_re > this->re_switch ) {
				pixel_val = this->_eval_function(this, x_ser, y_ser, r_ser, true);
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

			image[i + j*model->width] = scale * pixel_val;
			x += half_xbin;
		}
		y += half_ybin;
	}

}

/**
 * Constructor with sane defaults
 */
SersicLikeProfile::SersicLikeProfile() :
	Profile(),
	xcen(0), ycen(0),
	mag(15), box(0),
	ang(0), axrat(1),
	rough(false), acc(0.1),
	re_switch(1), resolution(9),
	max_recursions(2), adjust(true),
	re_max(0)
{
	// no-op
}

} /* namespace profit */
