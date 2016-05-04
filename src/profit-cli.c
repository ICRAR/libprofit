/**
 * Command-line utility to use profit
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "profit.h"
#include "sersic.h"
#include "sky.h"
#include "version.h"

char **_parse_profile_value(char *token) {
	char **key_and_val;
	char *equals = strchr(token, '=');
	if( !equals ) {
		fprintf(stderr, "Parameter %s doesn't give a value\n", token);
		return NULL;
	}
	if( strlen(equals) == 1 ) {
		fprintf(stderr, "Parameter %s gives an empty value\n", token);
		return NULL;
	}

	key_and_val = (char **)malloc(sizeof(char *) * 2);
	key_and_val[0] = strndup(token, equals - token);
	key_and_val[1] = strndup(equals + 1, strlen(token) - (token - equals));
	return key_and_val;
}

#define _READ_DOUBLE_OR_FAIL(key, val, name, len, dst) \
	do { \
		if ( !strncmp(key, name, len) ) { \
			tmp = strtod(val, &endptr); \
			if( tmp == 0 && endptr == val ) { \
				fprintf(stderr, "Invalid double value for %s: %s\n", key, val); \
				return -1;\
			} \
			dst = tmp;\
			return 1; \
		} \
	} while(0);

short _keyval_to_sersic(profit_profile *p, char *key, char *val) {

	char *endptr;
	double tmp;

	profit_sersic_profile *s = (profit_sersic_profile *)p;
	_READ_DOUBLE_OR_FAIL(key, val, "xcen",  4, s->xcen);
	_READ_DOUBLE_OR_FAIL(key, val, "ycen",  4, s->ycen);
	_READ_DOUBLE_OR_FAIL(key, val, "mag",   3, s->mag);
	_READ_DOUBLE_OR_FAIL(key, val, "re",    2, s->re);
	_READ_DOUBLE_OR_FAIL(key, val, "nser",  4, s->nser);
	_READ_DOUBLE_OR_FAIL(key, val, "box",   3, s->box);
	_READ_DOUBLE_OR_FAIL(key, val, "ang",   3, s->ang);
	_READ_DOUBLE_OR_FAIL(key, val, "axrat", 5, s->axrat);
	return 0;
}

short _keyval_to_sky(profit_profile *p, char *key, char *val) {
	char *endptr;
	double tmp;
	profit_sky_profile *s = (profit_sky_profile *)p;
	_READ_DOUBLE_OR_FAIL(key, val, "bg",  2, s->bg);
	return 0;
}

profit_profile *desc_to_profile(
	char *description,
	const char* name,
	unsigned short allow_empty_profile,
	short (keyval_to_param)(profit_profile *, char *, char *)
) {

	char *tok;
	char **key_and_val;
	short assigned;
	profit_profile *p;

	if( !description && !allow_empty_profile ) {
		fprintf(stderr, "Empty %s profile description\n", name);
		return NULL;
	}

	p = profit_get_profile(name);
	if( !description ) {
		return p;
	}

	while( (tok = strtok(description, ":")) ) {

		key_and_val = _parse_profile_value(tok);
		if( !key_and_val ) {
			fprintf(stderr, "Invalid token found in %s profile description: %s\n", name, description);
			return NULL;
		}

		assigned = keyval_to_param(p, key_and_val[0], key_and_val[1]);
		if( assigned == -1 ) {
			return NULL;
		}
		if( !assigned ) {
			fprintf(stderr, "Ignoring unknown %s profile parameter: %s\n", name, key_and_val[0]);
		}
		free(key_and_val);

		/* Otherwise we'll always start strtok from the beginning */
		description = NULL;
	}
	return p;
}

profit_profile *parse_profile(char *description) {

	/* The description might be only a name */
	char *subdesc = NULL;
	size_t name_end = strlen(description);
	char *colon = strchr(description, ':');
	if( colon ) {
		name_end = colon - description;
		subdesc = colon + 1;
	}

	if( !strncmp(description, "sersic", name_end) ) {
		return desc_to_profile(subdesc, "sersic", 0, &_keyval_to_sersic);
	}
	else if( !strncmp(description, "sky", name_end) ) {
		return desc_to_profile(subdesc, "sky", 0, &_keyval_to_sky);
	}

	fprintf(stderr, "Unknown profile name in profile description: %s\n", description);
	return NULL;

}

void usage(FILE *file, char *argv[]) {
	fprintf(file,"\n%s: utility program to generate an image out of a model and a set of profiles\n\n", argv[0]);
	fprintf(file,"This program is licensed under the GPLv3 license.\n\n");
	fprintf(file,"Usage: %s [options] -p <spec> [-p <spec> ...]\n\n",argv[0]);
	fprintf(file,"Options:\n");
	fprintf(file,"  -w        Image width. Defaults to 100\n");
	fprintf(file,"  -H        Image height. Defaults to 100\n\n");
	fprintf(file,"  -m        Zero magnitude. Defaults to 0.\n\n");
	fprintf(file,"  -h,-?     Show this help and exit\n");
	fprintf(file,"  -v        Show the program version and exit\n\n");
	fprintf(file,"Profiles should be specified as follows:\n\n");
	fprintf(file,"-p name:param1=val1:param2=val2:...\n\n");
}

typedef enum _output_type {
	none = 0,
	binary = 1,
	text = 2,
	csv = 3
} output_t;

int main(int argc, char *argv[]) {

	int opt;
	unsigned int width = 100, height = 100;
	double magzero = 0;
	unsigned int n_profiles = 0, i, j;
	char *endptr;
	output_t output = none;
	profit_profile *profile;
	profit_profile **profiles;

	/* Prepare for these many profiles */
	for(i=0; i!= argc; i++) {
		if( !strncmp("-p", argv[i], 2) ) {
			n_profiles++;
		}
	}
	if( !n_profiles ) {
		usage(stderr, argv);
		return 1;
	}

	profiles = (profit_profile **)malloc(sizeof(profit_profile *) * n_profiles);
	n_profiles = 0;

	while( (opt = getopt(argc, argv, "h?vp:w:H:m:tb")) != -1 ) {
		switch(opt) {

			case 'h':
			case '?':
				usage(stdout, argv);
				return 0;

			case 'v':
				printf("libprofit version %s\n", PROFIT_VERSION);
				return 0;

			case 'p':
				profile = parse_profile(optarg);
				if( profile == NULL ) {
					return 1;
				}
				profiles[n_profiles++] = profile;
				break;

			case 'w':
				width = atoi(optarg);
				break;

			case 'H':
				height = atoi(optarg);
				break;

			case 'm':
				magzero = strtod(optarg, &endptr);
				if( magzero == 0 && endptr == optarg ) {
					fprintf(stderr, "Invalid magzero value: %s\n", optarg);
					return 1;
				}
				break;

			case 't':
				if( output != none ) {
					fprintf(stderr, "-t and -b cannot be used together\n");
					return 1;
				}
				output = text;
				break;

			case 'b':
				if( output != none ) {
					fprintf(stderr, "-b and -t cannot be used together\n");
					return 1;
				}
				output = binary;
				break;

			default:
				usage(stderr, argv);
				return 1;

		}
	}

	/* We default to text output */
	if( output == none ) {
		output = text;
	}

	profit_model *m = (profit_model *)malloc(sizeof(profit_model));
	m->error      = NULL;
	m->width      = width;
	m->height     = height;
	m->res_x      = width;
	m->res_y      = height;
	m->magzero    = magzero;
	m->n_profiles = n_profiles;
	m->profiles   = profiles;

	/* Go, go, go */
	profit_make_model(m);

	/* Check for any errors */
	if( m->error ) {
		fprintf(stderr, "Error while calculating model: %s\n", m->error);
		return 1;
	}
	for(i=0; i!=n_profiles; i++) {
		if( m->profiles[i]->error ) {
			fprintf(stderr, "Error while calculating model: %s\n", m->profiles[i]->error);
			profit_cleanup(m);
			return 1;
		}
	}

	switch(output) {
		case binary:
			fwrite(m->image, sizeof(double), m->width * m->height, stdout);
			break;
		case text:
			for(j=0; j!=m->height; j++) {
				for(i=0; i!=m->width; i++) {
					printf("%g ", m->image[j*m->width + i]);
				}
				printf("\n");
			}
	}

	profit_cleanup(m);

}
