/**
 * Header file for moffat profile implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Aaron Robotham
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
#ifndef _MOFFAT_H_
#define _MOFFAT_H_

#include "sersic_like.h"

namespace profit
{

class MoffatProfile : public SersicLikeProfile {

protected:
	double get_lumtot(double r_box);
	double get_re();
	double adjust_re_switch();
	double adjust_re_max();
	eval_function_t get_evaluation_function();

public:

	MoffatProfile();

	/* General parameters */
	double fwhm;
	double con;

};

} /* namespace profit */

#endif /* _MOFFAT_H_ */
