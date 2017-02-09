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

typedef struct _d_point {
	double x;
	double y;
} d_point_t;

typedef struct _d_subsampling_info {
	d_point_t point;
	double xbin;
	double ybin;
	unsigned int resolution;
	unsigned int max_recursion;
} d_subsampling_info;

inline double d_evaluate_sersic(double x, double y, double box, double nser, double rscale, double bn) {
	private double r = pow(pow(fabs(x), 2+box) + pow(fabs(y), 2+box), 1/(2+box));
	private double r_factor = pow(r/rscale, 1/nser);
	return exp(-bn * (r_factor - 1));
}

inline void d_image_to_profile_coordiates(double x, double y, double *x_prof, double *y_prof, double xcen, double ycen, double cos_ang, double sin_ang, double axrat) {
	x -= xcen;
	y -= ycen;
	*x_prof =   x * cos_ang + y * sin_ang;
	*y_prof = (-x * sin_ang + y * cos_ang)/axrat;
}

kernel void sersic_double(
	global double *image,
	global d_point_t *to_subsample,
	int width, int height,
	int rough,
	double pixel_scale,
	double scale_x, double scale_y,
	double xcen, double ycen,
	double cos_ang, double sin_ang, double axrat,
	double rscale, double rscale_switch, double rscale_max,
	double box, double nser, double bn) {

	private int i = get_global_id(0);
	private double x = (i%width + 0.5)*scale_x;
	private double y = (i/width + 0.5)*scale_y;

	// image to profile coordinate conversion
	private double x_prof, y_prof;
	d_image_to_profile_coordiates(x, y, &x_prof, &y_prof, xcen, ycen, cos_ang, sin_ang, axrat);

	private double r_prof = sqrt(x_prof*x_prof + y_prof*y_prof);
	private double pixel_val;

	if( rscale_max > 0 && (r_prof/rscale) > rscale_max ) {
		pixel_val = 0.;
#if __OPENCL_C_VERSION__ <= 120
		to_subsample[i].x = -1;
#endif /* __OPENCL_C_VERSION__ */
	}
	else if( rough || (r_prof/rscale) > rscale_switch ) {
		pixel_val = d_evaluate_sersic(x_prof, y_prof, box, nser, rscale, bn);
#if __OPENCL_C_VERSION__ <= 120
		to_subsample[i].x = -1;
#endif /* __OPENCL_C_VERSION__ */
	}
	else {
		// subsample
		to_subsample[i].x = x;
		to_subsample[i].y = y;
	}

	image[i] = pixel_scale * pixel_val;

}

kernel void sersic_subsample_double(
	global double *image,
   global d_subsampling_info *all_info,
	double acc,
	double pixel_scale,
	double xcen, double ycen,
	double cos_ang, double sin_ang, double axrat,
	double rscale, double rscale_switch, double rscale_max,
	double box, double nser, double bn) {

	private int i = get_global_id(0);
	private d_subsampling_info info = all_info[i];
	private double x = info.point.x;
	private double y = info.point.y;

	// image to profile coordinate conversion
	// including delta_y_prof to test accuracy
	private double x_prof, y_prof;
	d_image_to_profile_coordiates(x, y, &x_prof, &y_prof, xcen, ycen, cos_ang, sin_ang, axrat);
	private double delta_y_prof = (-info.xbin * sin_ang + info.ybin * cos_ang)/axrat;

	private double val, testval;

	val = d_evaluate_sersic(x_prof, y_prof, box, nser, rscale, bn);
	testval = d_evaluate_sersic(x_prof, fabs(y_prof) + fabs(delta_y_prof), box, nser, rscale, bn);

	// no need for subsampling
	if( fabs(testval/val - 1.0) <= acc ) {
		all_info[i].point.x = -1.;
		all_info[i].point.y = -1.;
	}
	// else we already have the correct coordinates for the next subsampling

	image[i] = pixel_scale * val;

}
)==="
