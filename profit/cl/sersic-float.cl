R"===(
/**
 * Single-precision Sersic profile OpenCL kernel implementation for libprofit
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

typedef struct _f_point {
	float x;
	float y;
} f_point_t;

typedef struct _f_subsampling_info {
	f_point_t point;
	float xbin;
	float ybin;
	unsigned int resolution;
	unsigned int max_recursion;
} f_subsampling_info;


inline float f_evaluate_sersic(float x, float y, float box, float nser, float rscale, float bn) {
	private float r = pow(pow(fabs(x), 2+box) + pow(fabs(y), 2+box), 1/(2+box));
	private float r_factor = pow(r/rscale, 1/nser);
	return exp(-bn * (r_factor - 1));
}

inline void f_image_to_profile_coordiates(float x, float y, float *x_prof, float *y_prof, float xcen, float ycen, float cos_ang, float sin_ang, float axrat) {
	x -= xcen;
	y -= ycen;
	*x_prof =   x * cos_ang + y * sin_ang;
	*y_prof = (-x * sin_ang + y * cos_ang)/axrat;
}

kernel void sersic_float(
	global float *image,
	global f_point_t *to_subsample,
	int width, int height,
	int rough,
	float pixel_scale,
	float scale_x, float scale_y,
	float xcen, float ycen,
	float cos_ang, float sin_ang, float axrat,
	float rscale, float rscale_switch, float rscale_max,
	float box, float nser, float bn) {

	private int i = get_global_id(0);
	private float x = (i%width + 0.5f)*scale_x;
	private float y = (i/width + 0.5f)*scale_y;

	// image to profile coordinate conversion
	private float x_prof, y_prof;
	f_image_to_profile_coordiates(x, y, &x_prof, &y_prof, xcen, ycen, cos_ang, sin_ang, axrat);

	private float r_prof = sqrt(x_prof*x_prof + y_prof*y_prof);
	private float pixel_val;

	if( rscale_max > 0 && (r_prof/rscale) > rscale_max ) {
		pixel_val = 0.f;
#if __OPENCL_C_VERSION__ <= 120
		to_subsample[i].x = -1;
#endif /* __OPENCL_C_VERSION__ */
	}
	else if( rough || (r_prof/rscale) > rscale_switch ) {
		pixel_val = f_evaluate_sersic(x_prof, y_prof, box, nser, rscale, bn);
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

kernel void sersic_subsample_float(
	global float *image,
   global f_subsampling_info *all_info,
	float acc,
	float pixel_scale,
	float xcen, float ycen,
	float cos_ang, float sin_ang, float axrat,
	float rscale, float rscale_switch, float rscale_max,
	float box, float nser, float bn) {

	private int i = get_global_id(0);
	private f_subsampling_info info = all_info[i];
	private float x = info.point.x;
	private float y = info.point.y;

	// image to profile coordinate conversion
	// including delta_y_prof to test accuracy
	private float x_prof, y_prof;
	f_image_to_profile_coordiates(x, y, &x_prof, &y_prof, xcen, ycen, cos_ang, sin_ang, axrat);
	private float delta_y_prof = (-info.xbin * sin_ang + info.ybin * cos_ang)/axrat;

	private float val, testval;

	val = f_evaluate_sersic(x_prof, y_prof, box, nser, rscale, bn);
	testval = f_evaluate_sersic(x_prof, fabs(y_prof) + fabs(delta_y_prof), box, nser, rscale, bn);

	// no need for subsampling
	if( fabs(testval/val - 1.0f) <= acc ) {
		all_info[i].point.x = -1.f;
		all_info[i].point.y = -1.f;
	}
	// else we already have the correct coordinates for the next subsampling

	image[i] = pixel_scale * val;

}

)==="
