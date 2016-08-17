/**
 * Ferrer profile implementation
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

#include "ferrer.h"
#include "utils.h"

using namespace std;

namespace profit
{

/*
 * The evaluation of the ferrer profile at ferrer coordinates (x,y).
 *
 * The ferrer profile has this form:
 *
 * (1-r_factor)^(a)
 *
 * where r_factor = (r/re)^(2-b)
 *              r = (x^{2+B} + y^{2+B})^{1/(2+B)}
 *              B = box parameter
 */


/*
 * The main ferrer evaluation function for a given X/Y coordinate
 */
static inline
double _ferrer_for_xy_r(SersicLikeProfile *sp,
                        double x, double y,
                        double r, bool reuse_r) {

	FerrerProfile *fp = static_cast<FerrerProfile *>(sp);
	double r_factor;
	if( reuse_r && fp->box == 0 ) {
		r_factor = r;
	}
	else if( fp->box == 0 ) {
		r_factor = sqrt(x*x + y*y);
	}
	else {
		double box = fp->box + 2.;
		r_factor = pow( pow(abs(x), box) + pow(abs(y), box), 1./box);
	}

	r_factor /= fp->rout;
	if( r_factor < 1 ) {
		return pow(1 - pow(r_factor, 2 - fp->b), fp->a);
	}
	else {
		return 0;
	}
}

eval_function_t FerrerProfile::get_evaluation_function() {
	return &_ferrer_for_xy_r;
}

double FerrerProfile::get_lumtot(double r_box) {
	double a = this->a;
	double b = this->b;
	double g_factor = gammafn(a+1) * gammafn((4-b)/(2-b)) / gammafn(a+2/(2-b)+1);
	return pow(rout, 2) * M_PI * g_factor * axrat/r_box;
}

double FerrerProfile::adjust_re_switch() {
	return this->rout;
}

double FerrerProfile::adjust_re_max() {
	return this->rout;
}

double FerrerProfile::get_re() {
	return this->rout;
}

/**
 * Constructor with sane defaults
 */
FerrerProfile::FerrerProfile() :
	SersicLikeProfile(),
	rout(3), a(1), b(1)
{
	// no-op
}

} /* namespace profit */
