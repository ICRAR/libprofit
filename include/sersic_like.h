/**
 * Header file for sersic-like profile implementations
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
#ifndef _SERSIC_LIKE_H_
#define _SERSIC_LIKE_H_

#include "profit.h"

namespace profit
{

class SersicLikeProfile;
typedef double (*eval_function_t)(SersicLikeProfile *, double, double, double, bool);

class SersicLikeProfile : public Profile {

protected:

	void initial_calculations();
	void subsampling_params(double x, double y, unsigned int &res, unsigned int &max_rec);
	double get_pixel_scale();
	double adjust_acc();
	virtual double get_lumtot(double r_box) = 0;
	virtual double get_re() = 0;
	virtual double adjust_re_switch() = 0;
	virtual double adjust_re_max() = 0;
	virtual eval_function_t get_evaluation_function() = 0;

	/* These are internally calculated at profile evaluation time */
	double _ie;
	double _cos_ang;
	double _sin_ang;
	eval_function_t _eval_function;

private:

	void _image_to_profile_coordinates(double x, double y, double &x_prof, double &y_prof);

	double subsample_pixel(double x0, double x1,
	                       double y0, double y1,
	                       unsigned int recur_level,
	                       unsigned int max_recursions,
	                       unsigned int resolution);

public:

	SersicLikeProfile();
	void validate();
	void evaluate(double *image);

	/* General parameters */
	double xcen;
	double ycen;
	double mag;
	double ang;
	double axrat;
	double box;

	/*
	 * Common "Re" concept, profiles provide it in different ways
	 * via get_re()
	 */
	double _re;

	/* Used to control the subsampling */
	bool rough;
	double acc;
	double re_switch;
	unsigned int resolution;
	unsigned int max_recursions;
	bool adjust;

	/* Used to avoid outer regions */
	double re_max;

};

} /* namespace profit */

#endif /* _SERSIC_LIKE_H_ */
