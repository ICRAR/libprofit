R"===(
/**
 * float 2D convolution OpenCL implementation for libprofit
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

__kernel void convolve_float(
	const __global float *src,
	const int src_w,
	const int src_h,
	__constant float *krn,
	const __private int krn_w,
	const __private int krn_h,
	__global float *output
)
{

	const int X = get_global_id(0);
	const int Y = get_global_id(1);
	const int W = get_global_size(0);
	const int H = get_global_size(1);
	const int half_krn_w = krn_w / 2;
	const int half_krn_h = krn_h / 2;

	/* No need to do these */
	if (X >= src_w || Y >= src_h) {
		return;
	}

	/* perform convolution */
	__constant float *krn_ptr = krn + krn_w * krn_h - 1;
	const __global float *src_ptr = src + X - half_krn_w  + (Y - half_krn_h) * W;
	float sum = 0.0f;
	for (int j = 0; j < krn_h; j++) {
		int r = Y + j - half_krn_h;
		for (int i = 0; i < krn_w; i++) {
			int c = X + i - half_krn_w;
			if (c >= 0 && c < src_w &&
			    r >= 0 && r < src_h) {
				sum += *src_ptr * *krn_ptr;
			}
			src_ptr++;
			krn_ptr--;
		}
		src_ptr += W - krn_w;
	}

	output[X  + Y * W] = sum;
}

)==="
