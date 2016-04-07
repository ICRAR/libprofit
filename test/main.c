/**
 * Example usage of libprofit
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2014
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "profit.h"
#include "sersic.h"
#include "sky.h"

#define MATHLIB_STANDALONE 1
#include <Rmath.h>
double _Rf_qgamma_wrapper(double a, double b, double c) {
	return qgamma(a, b, c, 1, 0);
}
double _Rf_gammafn_wrapper(double a) {
	return gammafn(a);
}
double _Rf_beta_wrapper(double a, double b) {
	return beta(a, b);
}

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
double _gsl_qgamma_wrapper(double a, double b, double c) {
	return gsl_cdf_gamma_Qinv(a, b, c);
}
double _gsl_gammafn_wrapper(double a) {
	return gsl_sf_gamma(a);
}
double _gsl_beta_wrapper(double a, double b) {
	return gsl_sf_beta(a, b);
}

int main(int argc, char *argv[]) {

	profit_profile *buldge_p = profit_get_profile("sersic");
	profit_profile *disk_p = profit_get_profile("sersic");
	profit_profile *sky_p = profit_get_profile("sky");

	profit_sersic_profile *buldge_sp = (profit_sersic_profile *)buldge_p;
	buldge_sp->xcen  = 50.374;
	buldge_sp->ycen  = 93.229;
	buldge_sp->mag   = 17.70657;
	buldge_sp->re    = 4.252094;
	buldge_sp->nser  = 1.9602;

	profit_sersic_profile *disk_sp = (profit_sersic_profile *)disk_p;
	disk_sp->xcen  = 50.374;
	disk_sp->ycen  = 93.229;
	disk_sp->mag   = 17.70657;
	disk_sp->re    = 8.504189;
	disk_sp->nser  = 1.0;
	disk_sp->ang   = 9.5087;
	disk_sp->axrat = 0.406;

	profit_sky_profile *sky_sp = (profit_sky_profile *)sky_p;
	sky_sp->bg = 1.4e-16;

	profit_model *m = profit_get_model(3, sky_p, buldge_p, disk_p);
	m->width   = 101;
	m->height  = 185;
	m->res_x   = 101;
	m->res_y   = 185;
	m->magzero = 0.0;

	char *lib = argv[1];
	if( !strncmp(lib, "R", 1) ) {
		buldge_sp->_qgamma  = &_Rf_qgamma_wrapper;
		buldge_sp->_gammafn = &_Rf_gammafn_wrapper;
		buldge_sp->_beta = &_Rf_beta_wrapper;
		disk_sp->_qgamma  = &_Rf_qgamma_wrapper;
		disk_sp->_gammafn = &_Rf_gammafn_wrapper;
		disk_sp->_beta = &_Rf_beta_wrapper;
		set_seed(time(NULL), time(NULL));;
	}
	else if( !strncmp(lib, "gsl", 3) ) {
		buldge_sp->_qgamma  = &_gsl_qgamma_wrapper;
		buldge_sp->_gammafn = &_gsl_gammafn_wrapper;
		buldge_sp->_beta = &_gsl_beta_wrapper;
		disk_sp->_qgamma  = &_gsl_qgamma_wrapper;
		disk_sp->_gammafn = &_gsl_gammafn_wrapper;
		disk_sp->_beta = &_gsl_beta_wrapper;
	}

	int n = atoi(argv[2]);
	for(unsigned int i=0; i!=n; i++) {
		if( profit_make_model(m) ) {
			fputs("Error while calculating stuff", stderr);
			exit(1);
		}
		if( i == n-1 && argc <= 3 ) {
			free(m->image);
		}
	}

	if( argc > 3 ) {
		fwrite(m->image, sizeof(double), m->width * m->height, stdout);
		free(m->image);
	}

	free(m->profiles);

}
