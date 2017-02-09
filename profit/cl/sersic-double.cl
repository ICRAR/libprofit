R"===(
/**
 * Double-precision Sersic profile OpenCL kernel implementation for libprofit
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2017
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

#if __OPENCL_C_VERSION__ < 120
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#endif

inline double evaluate_sersic_double(double x, double y, double box, double nser, double rscale, double bn) {
	private double r = pow(pow(fabs(x), 2+box) + pow(fabs(y), 2+box), 1/(2+box));
	private double r_factor = pow(r/rscale, 1/nser);
	return exp(-bn * (r_factor - 1));
}

kernel void sersic_double(
	global double *image,
	int width, int height,
	double xcen, double ycen,
	double cos_ang, double sin_ang, double axrat,
	double rscale, double rscale_switch, double rscale_max, int rough,
	double box, double scale, double nser, double bn) {

	private int i = get_global_id(0);
	private double x = i % width + 0.5;
	private double y = i / width + 0.5;

	// image to profile coordinate conversion
	x -= xcen;
	y -= ycen;
	private double x_prof =  x * cos_ang + y * sin_ang;
	private double y_prof = -x * sin_ang + y * cos_ang;
	y_prof /= axrat;

	private double r_prof = sqrt(x_prof*x_prof + y_prof*y_prof);
	private double pixel_val;

	if( rscale_max > 0 && (r_prof/rscale) > rscale_max ) {
		pixel_val = 0.;
	}
	else if( rough || (r_prof/rscale) > rscale_switch ) {
		pixel_val = evaluate_sersic_double(x_prof, y_prof, box, nser, rscale, bn);
	}
	else {
		// subsample
		pixel_val = evaluate_sersic_double(x_prof, y_prof, box, nser, rscale, bn);
	}

	image[i] = scale * pixel_val;

}
)==="
