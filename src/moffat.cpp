/**
 * Moffat profile implementation
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

#include "moffat.h"
#include "utils.h"

using namespace std;

namespace profit
{

/*
 * The evaluation of the moffat profile at moffat coordinates (x,y).
 *
 * The moffat profile has this form:
 *
 * (1+r_factor^2)^(-c)
 *
 * where r_factor = (r/re)
 *              r = (x^{2+b} + y^{2+b})^{1/(2+b)}
 *              b = box parameter
 *
 * Reducing:
 *  r_factor = ((x/re)^{2+b} + (y/re)^{2+b})^{1/(2+b)}
 */

/*
 * The main moffat evaluation function for a given X/Y coordinate
 */
double _moffat_for_xy_r(SersicLikeProfile *sp,
                        double x, double y,
                        double r, bool reuse_r) {

	MoffatProfile *mp = static_cast<MoffatProfile *>(sp);
	double r_factor = (mp->box == 0) ?
	                  sqrt(x*x + y*y)/mp->_re :
	                  (pow(pow(abs(x),2.+mp->box)+pow(abs(y),2.+mp->box),1./(2.+mp->box)) ) / mp->_re;

	return 1/(pow(1+pow(r_factor,2), mp->con));
}

eval_function_t MoffatProfile::get_evaluation_function() {
	return &_moffat_for_xy_r;
}

double MoffatProfile::get_lumtot(double r_box) {
	double fwhm = this->fwhm;
	double con = this->con;
	return pow(this->_re, 2) * M_PI * axrat/(con-1)/r_box;
}

double MoffatProfile::get_re() {
	return fwhm/(2*sqrt(pow(2,(1/con))-1));
}

double MoffatProfile::adjust_re_switch() {
	double re_switch = this->fwhm*4;
	re_switch = max(min(re_switch, 20.), 2.);
	return re_switch / this->_re;
}

double MoffatProfile::adjust_re_max() {
	return this->fwhm*8;
}

/**
 * Constructor with sane defaults
 */
MoffatProfile::MoffatProfile() :
	SersicLikeProfile(),
	fwhm(3), con(2)
{
	// no-op
}

} /* namespace profit */
