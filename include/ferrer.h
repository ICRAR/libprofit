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

#include "profit.h"

namespace profit
{

class FerrerProfile : public Profile {

public:

	FerrerProfile();

	void validate();
	void evaluate(double *image);

	/* General parameters */
	double xcen;
	double ycen;
	double mag;
	double rout;
	double a;
	double b;
	double ang;
	double axrat;
	double box;

	/* Used to control the subsampling */
	bool rough;
	double acc;
	double re_switch;
	unsigned int resolution;
	unsigned int max_recursions;
	bool adjust;

	/* Used to avoid outer regions */
	double re_max;

	/* These are internally calculated profile init */
	double _ie;
	double _cos_ang;
	double _sin_ang;

};

} /* namespace profit */

#endif /* _FERRER_H_ */
