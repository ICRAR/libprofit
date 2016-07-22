/**
 * libprofit main entry-point routines
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

#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "convolve.h"
#include "profit.h"
#include "psf.h"
#include "sersic.h"
#include "sky.h"
#include "utils.h"

namespace profit {

struct _profit_profile_index {
	char *name;
	Profile *(* create)(void);
};

static
struct _profit_profile_index _all_profiles[] = {
	{"sky",    profit_create_sky},
	{"sersic", profit_create_sersic},
	{"psf",    profit_create_psf},
	{NULL, NULL} // Sentinel
};

Model::Model() :
	width(0), height(0),
	res_x(0), res_y(0),
	magzero(0),
	psf(NULL), psf_width(0), psf_height(0),
	calcmask(NULL), image(NULL),
	profiles(),
	error(NULL)
{
	// no-op
}

Profile* Model::add_profile(const char * profile_name) {

	struct _profit_profile_index *p = _all_profiles;
	while(1) {
		if( p->name == NULL ) {
			break;
		}
		if( !strcmp(profile_name, p->name) ) {
			Profile *profile = p->create();
			profile->model = this;
			profile->error = NULL;
			profile->name = profile_name;
			profile->convolve = false;
			this->profiles.push_back(profile);
			return profile;
		}
		p++;
	}

	return NULL;
}

void Model::evaluate() {

	unsigned int p;

	/* Check limits */
	if( !this->width ) {
		this->error = strdup("Model's width is 0");
		return;
	}
	else if( !this->height ) {
		this->error = strdup("Model's height is 0");
		return;
	}
	else if( !this->res_x ) {
		this->error = strdup("Model's res_x is 0");
		return;
	}
	else if( !this->res_y ) {
		this->error = strdup("Model's res_y is 0");
		return;
	}

	/*
	 * If at least one profile is requesting convolving we require
	 * a valid psf.
	 */
	std::vector<Profile *>::iterator pit;
	for(p=0, pit=this->profiles.begin(); pit!=this->profiles.end(); pit++, p++) {
		Profile *profile = *pit;
		if( profile->convolve ) {
			if( !this->psf ) {
				const char *msg = "Profile %s requires convolution but no psf was provided";
				this->error = (char *)malloc(strlen(msg) - 1 + strlen(profile->name));
				sprintf(this->error, msg, profile->name);
				return;
			}
			if( !this->psf_width ) {
				this->error = strdup("Model's psf width is 0");
				return;
			}
			if( !this->psf_height ) {
				this->error = strdup("Model's psf height is 0");
				return;
			}
			break;
		}
	}

	this->xbin = this->width/(double)this->res_x;
	this->ybin = this->height/(double)this->res_y;
	this->image = (double *)calloc(this->width * this->height, sizeof(double));
	if( !this->image ) {
		char *msg = "Cannot allocate memory for image with w=%u, h=%u";
		this->error = (char *)malloc( strlen(msg) - 4 + 20 ); /* 32bits unsigned max is 4294967295 (10 digits) */
		sprintf(this->error, msg, this->width, this->height);
		return;
	}

	/*
	 * Validate all profiles.
	 * Each profile can fail during validation in which case we don't proceed any further
	 */
	for(p=0, pit=this->profiles.begin(); pit!=this->profiles.end(); pit++, p++) {
		Profile *profile = *pit;
		profile->validate();
		if( profile->error ) {
			return;
		}
	}

	/*
	 * Generate a separate image for each profile.
	 *
	 * We optionally use OpenMP to parallelize this per-profile image
	 * generation. Depending on how many there are we might get a speed up, so
	 * probably we should study what is the best way to go here (e.g.,
	 * parallelize only if we have more than 2 or 3 profiles)
	 */
	double **profile_images = (double **)malloc(sizeof(double *) * this->profiles.size());
#if _OPENMP
	#pragma omp parallel for private(p)
#endif
	for(p=0, pit=this->profiles.begin(); pit!=this->profiles.end(); pit++, p++) {
		Profile *profile = *pit;
		profile_images[p] = (double *)calloc(this->width * this->height, sizeof(double));
		profile->evaluate(profile_images[p]);
	}

	/*
	 * Sum up all results
	 *
	 * We first sum up all images that need convolving, we convolve them
	 * and after that we add up the remaining images.
	 */
	bool convolve = false;
	for(p=0, pit=this->profiles.begin(); pit!=this->profiles.end(); pit++, p++) {
		Profile *profile = *pit;
		if( profile->convolve ) {
			convolve = true;
			profit_add_images(this->image, profile_images[p], this->width, this->height);
		}
	}
	if( convolve ) {
		size_t psf_size = sizeof(double) * this->psf_width * this->psf_height;
		double *psf = (double *)malloc(psf_size);
		memcpy(psf, this->psf, psf_size);
		profit_normalize(psf, this->psf_width, this->psf_height);
		profit_convolve(this->image, this->width, this->height, psf, this->psf_width, this->psf_height, this->calcmask, true);
		free(psf);
	}
	for(p=0, pit=this->profiles.begin(); pit!=this->profiles.end(); pit++, p++) {
		Profile *profile = *pit;
		if( !profile->convolve ) {
			profit_add_images(this->image, profile_images[p], this->width, this->height);
		}
		free(profile_images[p]);
	}
	free(profile_images);

	/* Done! Good job :-) */
}

char *Model::get_error() {

	if( this->error ) {
		return this->error;
	}

	std::vector<Profile *>::iterator p;
	for(p=this->profiles.begin(); p!=this->profiles.end(); p++) {
		Profile *profile = *p;
		if( profile->error ) {
			return profile->error;
		}
	}
	return NULL;
}

Model::~Model() {

	std::vector<Profile *>::iterator p;
	for(p=this->profiles.begin(); p!=this->profiles.end(); p++) {
		Profile *profile = *p;
		free(profile->error);
		delete profile;
	}
	free(this->error);
	free(this->image);
	free(this->psf);
	free(this->calcmask);
}

} /* namespace profit */
