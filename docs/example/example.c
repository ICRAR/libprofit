/* copyright statement, etc */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "example.h"

void profit_init_example(profit_profile *profile, profit_model *model) {

	profit_example_profile *e = (profit_example_profile *)profile;

	if ( e->param1 < 0 ) {
		profile->error = strdup("param1 is negative");
		return;
	}
	if ( e->param1 < 0 ) {
		profile->error = strdup("param2 is negative");
		return;
	}
	if ( e->param3 < 0 ) {
		profile->error = strdup("param3 is negative");
		return;
	}

	/*
	if ( model->width < 20 || model->height < 20 ) {
		profile->error = strdup("can't apply example profile to images less than 20x20");
		return;
	}
	*/
}

void profit_make_example(profit_profile *profile, profit_model *model, double *image) {

	double x, y;
	unsigned int i, j;
	double half_xbin = model->xbin/2.;
	double half_ybin = model->ybin/2.;

	profit_example_profile *e = (profit_example_profile *)profile;

	x = 0;
	for (i=0; i < model->width; i++) {
		x += half_xbin;

		y = 0;
		for (j=0; j < model->height; j++) {
			y += half_ybin;

			double val = fabs( (e->param1 - e->param2) * e->param3 * (x - y) );
			image[i + j*model->width] = val;

			y += half_ybin;
		}
		x += half_xbin;
	}
}

profit_profile *profit_create_example() {

	profit_example_profile *e = (profit_example_profile*) malloc(sizeof(profit_example_profile));
	e->profile.init_profile = &profit_init_example;
	e->profile.make_profile = &profit_make_example;

	e->param1 = 1.;
	e->param2 = 2.;
	e->param3 = 3;

	return (profit_profile *)e;
}
