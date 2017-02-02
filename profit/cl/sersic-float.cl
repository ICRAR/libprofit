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

inline float evaluate_sersic_float(float x, float y, float box, float nser, float rscale, float bn) {
	private float r = pow(pow(fabs(x), 2+box) + pow(fabs(y), 2+box), 1/(2+box));
	private float r_factor = pow(r/rscale, 1/nser);
	return exp(-bn * (r_factor - 1));
}

kernel void sersic_float(global float *image,
	int width, int height,
	float xcen, float ycen,
	float cos_ang, float sin_ang, float axrat,
	float rscale, float rscale_switch, float rscale_max, int rough,
	float box, float scale, float nser, float bn) {

	private int i = get_global_id(0);
	private float x = i % width + 0.5f;
	private float y = i / width + 0.5f;

	// image to profile coordinate conversion
	x -= xcen;
	y -= ycen;
	private float x_prof =  x * cos_ang + y * sin_ang;
	private float y_prof = -x * sin_ang + y * cos_ang;
	y_prof /= axrat;

	private float r_prof = sqrt(x_prof*x_prof + y_prof*y_prof);
	private float pixel_val;

	if( rscale_max > 0 && (r_prof/rscale) > rscale_max ) {
		pixel_val = 0.f;
	}
	else if( rough || (r_prof/rscale) > rscale_switch ) {
		pixel_val = evaluate_sersic_float(x, y, box, nser, rscale, bn);
	}
	else {
		// subsample
		pixel_val = evaluate_sersic_float(x, y, box, nser, rscale, bn);
	}

	image[i] = scale * pixel_val;

}
)==="
