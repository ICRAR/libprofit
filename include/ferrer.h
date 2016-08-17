/**
 * Header file for ferrer profile implementation
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
#ifndef _FERRER_H_
#define _FERRER_H_

#include "sersic_like.h"

namespace profit
{

class FerrerProfile : public SersicLikeProfile {

protected:
	double get_lumtot(double r_box);
	double get_rscale();
	double adjust_rscale_switch();
	double adjust_rscale_max();
	double adjust_acc();
	eval_function_t get_evaluation_function();

public:

	FerrerProfile();


	/* General parameters */
	double rout;
	double a;
	double b;

};

} /* namespace profit */

#endif /* _FERRER_H_ */
