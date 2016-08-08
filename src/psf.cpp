/**
 * PSF profile implementation
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Aaron Robotham, Rodrigo Tobar
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

#include <algorithm>
#include <cmath>

#include "psf.h"
#include "utils.h"

using namespace std;

namespace profit
{

void PsfProfile::validate()  {

	if( !this->model->psf ) {
		throw invalid_parameter("No psf present in the model, cannot produce a psf profile");
	}
	this->scale = pow(10, -0.4*(this->mag - this->model->magzero));

}

static inline
void psf_apply(PsfProfile *psf, Model *model, double *image,
               double *psf_img, unsigned int psf_w, unsigned int psf_h,
               int target_x, int target_y) {

	unsigned int i, j, img_x, img_y;

	for(j=0; j!=psf_h; j++) {

		/* Don't draw outside the boundaries of the full image */
		if( (int)j+target_y < 0 ) {
			continue;
		}
		img_y = j + (unsigned int)target_y;
		if( img_y >= model->height ) {
			break;
		}

		for(i=0; i!=psf_w; i++) {

			/* Don't draw outside the boundaries of the full image */
			if( (int)i+target_x < 0 ) {
				continue;
			}
			img_x = i + (unsigned int)target_x;
			if( img_x >= model->width ) {
				break;
			}

			image[img_x + img_y*model->width] = psf_img[i + j*psf_w];
		}
	}

}

double * regrid_and_normalize(PsfProfile *psf_profile, Model *model,
                              unsigned int &psf_width,
                              unsigned int &psf_height) {

	/* if pixels are the same size simply return the psf as is */
	if( model->scale_x == model->psf_scale_x &&
	    model->scale_y == model->psf_scale_y ) {
		psf_width  = model->psf_width;
		psf_height = model->psf_height;
		double *psf = new double[psf_width * psf_height];
		copy(model->psf, model->psf + (psf_width * psf_height), psf);
		normalize(psf, psf_width, psf_height);
		return psf;
	}

	/* Otherwise re-grid the psf in model-pixel units */
	unsigned int w = psf_width = (unsigned int)ceil(model->psf_width  / model->scale_x);
	unsigned int h = psf_height = (unsigned int)ceil(model->psf_height / model->scale_y);
	double *regridded_psf = new double[w*h];

	for(unsigned j=0; j!=h; j++) {
		for(unsigned i=0; i!=w; i++) {

			// each pixel borrows a bit from the original image
			double value = 1;
			regridded_psf[i + j*w] = value;
		}
	}

	normalize(regridded_psf, w, h);
	return regridded_psf;
}

void PsfProfile::evaluate(double *image) {

	unsigned int i, j;

	/*
	 * First of all, translate the PSF from its own
	 * pixel sizes to the image's pixel sizes; otherwise the two
	 * are not comparable.
	 */
	unsigned int psf_width, psf_height;
	double *psf = regrid_and_normalize(this, model, psf_width, psf_height);

	/*
	 * The PSF is not simply put "as is" in the nearest position of the desired
	 * center. Its values are interpolated instead to take into account the
	 * sub-pixel distance needed to reach the center.
	 *
	 * We avoid this extra calculation only if the target position of the psf's
	 * origin is located exactly on a pixel crossing, and if the psf's
	 * dimensions are both even. This would ensure that each pixel of the PSF
	 * corresponds exactly to one pixel on the target image, allowing us to have
	 * a direct copy of values.
	 */
	double psf_origin_x = (this->xcen - psf_width/2.) / model->scale_x;
	double psf_origin_y = (this->ycen - psf_height/2.)/ model->scale_y;
	if( (psf_width % 2 == 0 && psf_height % 2 == 0) && \
	    (floor(psf_origin_x) == psf_origin_x || ceil(psf_origin_x) == psf_origin_x) && \
	    (floor(psf_origin_y) == psf_origin_y || ceil(psf_origin_y) == psf_origin_y) ) {

		psf_apply(this, model, image,
		          psf, psf_width, psf_height,
		          (int)psf_origin_x, (int)psf_origin_y);

		return;
	}

	/*
	 * Each image pixel is now divided into four regions because of its
	 * intersection with the PSF pixels:
	 *
	 *             |
	 *    ---------|------
	 *    |        |      |
	 *    |   a3   |  a4  |  yd2
	 *  ___________|________
	 *    |        |      |
	 *    |        |      |
	 *    |   a1   |  a2  |  yd1
	 *    |        |      |
	 *    ---------|-------
	 *        xd1  |  xd2
	 *
	 * We average the four areas to obtain the value of the pixel. The areas are
	 * all the same on each image pixel so we calculate them once.
	 */

	double xd1 = psf_origin_x - floor(psf_origin_x);
	double xd2 = 1 - xd1;
	double yd1 = psf_origin_y - floor(psf_origin_y);
	double yd2 = 1 - yd1;
	double a1 = xd1 * yd1;
	double a2 = xd2 * yd1;
	double a3 = xd1 * yd2;
	double a4 = xd2 * yd2;

	unsigned int new_psf_w = psf_width + 1;
	unsigned int new_psf_h = psf_height + 1;
	double *new_psf = new double[new_psf_w * new_psf_h];

	for(j=0; j!=new_psf_h; j++) {
		for(i=0; i!=new_psf_w; i++) {

			/* The borders of the target image area use less psf pixels */
			double psf_val = 0;
			if( i != 0 && j != 0 ) {
				psf_val += psf[i-1 + (j-1)*psf_width] * a1;
			}
			if( i != 0 && j != (new_psf_h - 1) ) {
				psf_val += psf[i-1 + j*psf_width] * a3;
			}
			if( i != (new_psf_w - 1) && j != 0 ) {
				psf_val += psf[i + (j-1)*psf_width] * a2;
			}
			if( i != (new_psf_w - 1) && j != (new_psf_h - 1) ) {
				psf_val += psf[i + j*psf_width] * a4;
			}

			new_psf[i + j*new_psf_w] = psf_val;
		}
	}

	psf_apply(this, model, image,
	          new_psf, new_psf_w, new_psf_h,
	          (int)floor(psf_origin_x), (int)floor(psf_origin_y));

	delete [] new_psf;

}

PsfProfile::PsfProfile() :
	Profile(),
	xcen(0),
	ycen(0),
	mag(0)
{
	// no-op
}

} /* namespace profit */
