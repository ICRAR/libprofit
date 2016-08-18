/**
 * Header file for analytical profile implementations
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
#ifndef _ANALYTIC_H_
#define _ANALYTIC_H_

#ifdef PROFIT_DEBUG
#include <map>
#endif


#include "profit.h"

namespace profit
{

class AnalyticProfile;
typedef double (*eval_function_t)(AnalyticProfile *, double, double, double, bool);

class AnalyticProfile : public Profile {

protected:

	void initial_calculations();
	void subsampling_params(double x, double y, unsigned int &res, unsigned int &max_rec);
	double get_pixel_scale();
	double adjust_acc();
	virtual double get_lumtot(double r_box) = 0;
	virtual double get_rscale() = 0;
	virtual double adjust_rscale_switch() = 0;
	virtual double adjust_rscale_max() = 0;
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

	AnalyticProfile();
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
	 * radius scale, profiles provide it in different ways
	 * via get_rscale()
	 */
	double rscale;

	/* Used to control the subsampling */
	bool rough;
	double acc;
	double rscale_switch;
	unsigned int resolution;
	unsigned int max_recursions;
	bool adjust;

	/* Used to avoid outer regions */
	double rscale_max;

#ifdef PROFIT_DEBUG
	/* record of how many subintegrations we've done */
	std::map<int,int> n_integrations;
#endif

};

} /* namespace profit */

#endif /* _ANALYTIC_H_ */
