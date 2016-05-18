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
#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "profit.h"
#include "psf.h"
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

	size_t equals_idx = (size_t)(equals - token);
	key_and_val = (char **)malloc(sizeof(char *) * 2);
	key_and_val[0] = strndup(token, equals_idx);
	key_and_val[1] = strndup(equals + 1, strlen(token) - equals_idx);
	return key_and_val;
}

#define _READ_DOUBLE_OR_FAIL(key, val, name, len, dst) \
	do { \
		char *endptr; \
		double tmp; \
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

#define _READ_BOOL_OR_FAIL(key, val, name, len, dst) \
	do { \
		char *endptr; \
		long int tmp; \
		if ( !strncmp(key, name, len) ) { \
			tmp = strtol(val, &endptr, 10); \
			if( tmp == 0 && endptr == val ) { \
				fprintf(stderr, "Invalid integer value for %s: %s\n", key, val); \
				return -1;\
			} \
			dst = (bool)tmp;\
			return 1; \
		} \
	} while(0);

short _keyval_to_sersic(profit_profile *p, char *key, char *val) {
	profit_sersic_profile *s = (profit_sersic_profile *)p;
	_READ_DOUBLE_OR_FAIL(key, val, "xcen",  4, s->xcen);
	_READ_DOUBLE_OR_FAIL(key, val, "ycen",  4, s->ycen);
	_READ_DOUBLE_OR_FAIL(key, val, "mag",   3, s->mag);
	_READ_DOUBLE_OR_FAIL(key, val, "re",    2, s->re);
	_READ_DOUBLE_OR_FAIL(key, val, "nser",  4, s->nser);
	_READ_DOUBLE_OR_FAIL(key, val, "box",   3, s->box);
	_READ_DOUBLE_OR_FAIL(key, val, "ang",   3, s->ang);
	_READ_DOUBLE_OR_FAIL(key, val, "axrat", 5, s->axrat);
	_READ_BOOL_OR_FAIL(key, val, "rough",    5, s->rough);
	_READ_BOOL_OR_FAIL(key, val, "convolve", 8, p->convolve);
	return 0;
}

short _keyval_to_sky(profit_profile *p, char *key, char *val) {
	profit_sky_profile *s = (profit_sky_profile *)p;
	_READ_DOUBLE_OR_FAIL(key, val, "bg",  2, s->bg);
	_READ_BOOL_OR_FAIL(key, val, "convolve", 8, p->convolve);
	return 0;
}

short _keyval_to_psf(profit_profile *p, char *key, char *val) {
	profit_psf_profile *s = (profit_psf_profile *)p;
	_READ_DOUBLE_OR_FAIL(key, val, "xcen",  4, s->xcen);
	_READ_DOUBLE_OR_FAIL(key, val, "ycen",  4, s->ycen);
	_READ_DOUBLE_OR_FAIL(key, val, "mag",   3, s->mag);
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

	p = profit_create_profile(name);
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
		free(key_and_val[0]);
		free(key_and_val[1]);
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
		name_end = (size_t)(colon - description);
		subdesc = colon + 1;
	}

	if( !strncmp(description, "sersic", name_end) ) {
		return desc_to_profile(subdesc, "sersic", 0, &_keyval_to_sersic);
	}
	else if( !strncmp(description, "sky", name_end) ) {
		return desc_to_profile(subdesc, "sky", 0, &_keyval_to_sky);
	}
	else if( !strncmp(description, "psf", name_end) ) {
		return desc_to_profile(subdesc, "psf", 0, &_keyval_to_psf);
	}

	fprintf(stderr, "Unknown profile name in profile description: %s\n", description);
	return NULL;

}

double *parse_psf(char *optarg, unsigned int *psf_width, unsigned int *psf_height) {

	char *tok, *values, *endptr;
	unsigned int size, i = 0;
	double *psf;

	/* format is w:h:val1,val2... */
	tok = strtok(optarg, ":");
	if( !tok ) {
		fprintf(stderr, "Missing psf's width\n");
		return NULL;
	}
	*psf_width = (unsigned int)strtoul(tok, &endptr, 10);
	if( tok == endptr ) {
		fprintf(stderr, "Invalid value for psf's width: %s\n", tok);
		return NULL;
	}

	tok = strtok(NULL, ":");
	if( !tok ) {
		fprintf(stderr, "Missing psf's height\n");
		return NULL;
	}
	*psf_height = (unsigned int)strtoul(tok, &endptr, 10);
	if( tok == endptr ) {
		fprintf(stderr, "Invalid value for psf's height: %s\n", tok);
		return NULL;
	}

	values = strtok(NULL, ":");
	if( !values ) {
		fprintf(stderr, "Missing psf's values\n");
		return NULL;
	}

	size = *psf_width * *psf_height;
	psf = (double *)malloc(sizeof(double) * size);
	while( (tok = strtok(values, ",")) ) {
		values = NULL;
		psf[i] = strtod(tok, &endptr);
		if( psf[i] == 0 && tok == endptr ) {
			fprintf(stderr, "Invalid floating-point value for psf: %s\n", tok);
			free(psf);
			return NULL;
		}
		i++;
	}

	if( i != size ) {
		fprintf(stderr, "Not enough values provided for PSF. Provided: %u, expected: %u\n", i, size);
		free(psf);
		return NULL;
	}

	return psf;
}

void usage(FILE *file, char *argv[]) {
	fprintf(file,"\n%s: utility program to generate an image out of a model and a set of profiles\n\n", argv[0]);
	fprintf(file,"This program is licensed under the GPLv3 license.\n\n");
	fprintf(file,"Usage: %s [options] -p <spec> [-p <spec> ...]\n\n",argv[0]);
	fprintf(file,"Options:\n");
	fprintf(file,"  -t        Output image as text values on stdout\n");
	fprintf(file,"  -b        Output image as binary content on stdout\n");
	fprintf(file,"  -f <file> Output image as fits file\n");
	fprintf(file,"  -i <n>    Output performance information after evaluating the model n times\n");
	fprintf(file,"  -w        Image width. Defaults to 100\n");
	fprintf(file,"  -H        Image height. Defaults to 100\n");
	fprintf(file,"  -m        Zero magnitude. Defaults to 0.\n");
	fprintf(file,"  -P        PSF function (specified as w:h:val1,val2...)\n");
	fprintf(file,"  -h,-?     Show this help and exit\n");
	fprintf(file,"  -v        Show the program version and exit\n\n");
	fprintf(file,"Profiles should be specified as follows:\n\n");
	fprintf(file,"-p name:param1=val1:param2=val2:...\n\n");
	fprintf(file,"The following profiles (and parameters) are currently accepted:\n\n");
	fprintf(file," * psf: xcen, ycen, mag\n");
	fprintf(file," * sky: bg\n");
	fprintf(file," * sersic: xcen, ycen, mag, re, nser, box, ang, axrat, rough\n\n");
}

static inline
bool is_little_endian() {
	volatile uint32_t i=0x01234567;
	return (*((uint8_t*)(&i))) == 0x67;
}

static inline
double to_bigendian(double v) {
	double r;
	char *vbytes = (char *)(&v);
	char *rbytes = (char *)(&r);
	rbytes[0] = vbytes[7];
	rbytes[1] = vbytes[6];
	rbytes[2] = vbytes[5];
	rbytes[3] = vbytes[4];
	rbytes[4] = vbytes[3];
	rbytes[5] = vbytes[2];
	rbytes[6] = vbytes[1];
	rbytes[7] = vbytes[0];
	return r;
}

int to_fits(profit_model *m, char *fits_output) {

	FILE *f;
	unsigned int i, j, pos, padding;
	char hdr[80];

	/* Append .fits if not in the name yet */
	if( strstr(fits_output, ".fits") != &fits_output[strlen(fits_output) - 5] ) {
		fits_output = (char *)realloc(fits_output, strlen(fits_output) + 5);
		strcat(fits_output, ".fits");
	}

	f = fopen(fits_output, "wb");
	free(fits_output);
	if( !f ) {
		profit_cleanup(m);
		return 1;
	}

	/* Standard headers*/
	fprintf(f, "%-80s", "SIMPLE  = T               / File conforms to FITS standard");
	fprintf(f, "%-80s", "BITPIX  = -64             / Bits per pixel");
	fprintf(f, "%-80s", "NAXIS   = 2               / Number of axes");
	sprintf(hdr, "NAXIS1  = %-10.0u                / Width", m->width);
	fprintf(f, "%-80s", hdr);
	sprintf(hdr, "NAXIS2  = %-10.0u                / Height", m->height);
	fprintf(f, "%-80s", hdr);
	fprintf(f, "%-80s", "CRPIX1  = 1");
	fprintf(f, "%-80s", "CRVAL1  = 1");
	fprintf(f, "%-80s", "CDELT1  = 1");
	fprintf(f, "%-80s", "CTYPE1  = ' '");
	fprintf(f, "%-80s", "CUNIT1  = ' '");
	fprintf(f, "%-80s", "CRPIX2  = 1");
	fprintf(f, "%-80s", "CRVAL2  = 1");
	fprintf(f, "%-80s", "CDELT2  = 1");
	fprintf(f, "%-80s", "CTYPE2  = ' '");
	fprintf(f, "%-80s", "CUNIT2  = ' '");
	fprintf(f, "%-80s", "END");

	pos = (unsigned int)ftell(f);
	padding = 36*80 - (pos%36*80);
	for(i=0; i<padding; i++) {
		fprintf(f, " ");
	}

	/* data has to be big-endian */
	if( is_little_endian() ) {
		double *big_endian_image = (double *)malloc(sizeof(double) * m->width * m->height);
		for(j=0; j!=m->height; j++) {
			for(i=0; i!=m->width; i++) {
				pos = i + j*m->width;
				big_endian_image[pos] = to_bigendian(m->image[pos]);
			}
		}

		/* Simply replace the model's image, nobody will use it but us */
		free(m->image);
		m->image = big_endian_image;
	}

	fwrite(m->image, sizeof(double), m->width * m->height, f);

	/* Pad with zeroes until we complete the current 36*80 block */
	padding = 36*80 - ((m->width * m->height) % (36*80));
	void *zeros = calloc(padding, sizeof(double));
	fwrite(zeros, sizeof(double), padding, f);
	free(zeros);
	fclose(f);

	return 0;
}

typedef enum _output_type {
	none = 0,
	binary = 1,
	text = 2,
	fits = 3,
	performance = 4
} output_t;

int main(int argc, char *argv[]) {

	int opt;
	unsigned int width = 100, height = 100, iterations = 1;
	long duration;
	double magzero = 0, *psf = NULL;
	unsigned int i, j, psf_width = 0, psf_height = 0;
	char *endptr, *error, *fits_output = NULL;
	output_t output = none;
	profit_profile *profile;
	profit_model *m = profit_create_model();

	while( (opt = getopt(argc, argv, "h?vP:p:w:H:m:tbf:i:")) != -1 ) {
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
					profit_cleanup(m);
					return 1;
				}
				profit_add_profile(m, profile);
				break;

			case 'P':
				psf = parse_psf(optarg, &psf_width, &psf_height);
				if( !psf ) {
					usage(stderr, argv);
					profit_cleanup(m);
					return 1;
				}
				break;

			case 'w':
				width = (unsigned int)atoi(optarg);
				break;

			case 'H':
				height = (unsigned int)atoi(optarg);
				break;

			case 'm':
				magzero = strtod(optarg, &endptr);
				if( magzero == 0 && endptr == optarg ) {
					fprintf(stderr, "Invalid magzero value: %s\n", optarg);
					profit_cleanup(m);
					return 1;
				}
				break;

			case 't':
				if( output != none ) {
					fprintf(stderr, "-t and -b cannot be used together\n");
					profit_cleanup(m);
					return 1;
				}
				output = text;
				break;

			case 'b':
				if( output != none ) {
					fprintf(stderr, "-b and -t cannot be used together\n");
					profit_cleanup(m);
					return 1;
				}
				output = binary;
				break;

			case 'f':
				free(fits_output);
				fits_output = strdup(optarg);
				output = fits;
				break;

			case 'i':
				iterations = (unsigned int)atoi(optarg);
				output = performance;
				break;

			default:
				usage(stderr, argv);
				profit_cleanup(m);
				return 1;

		}
	}

	/* No profiles given */
	if( !m->n_profiles ) {
		usage(stderr, argv);
		return 1;
	}

	/* We default to text output */
	if( output == none ) {
		output = text;
	}

	m->width      = width;
	m->height     = height;
	m->res_x      = width;
	m->res_y      = height;
	m->magzero    = magzero;
	m->psf        = psf;
	m->psf_width  = psf_width;
	m->psf_height = psf_height;

	/* Go, go, go */
	profit_eval_model(m);

	if( output == performance ) {

		/* This means that we evaluated the model once, but who cares */
		struct timespec start, end;
		clock_gettime(CLOCK_MONOTONIC, &start);
		for(i=0; i!=iterations; i++) {
			free(m->image);
			free(m->error);
			profit_eval_model(m);
		}
		clock_gettime(CLOCK_MONOTONIC, &end);
		duration = (end.tv_sec - start.tv_sec)*1000000 + (end.tv_nsec - start.tv_nsec)/1000;

	}

	/* Check for any errors */
	error = profit_get_error(m);
	if( error ) {
		fprintf(stderr, "Error while calculating model: %s\n", error);
		profit_cleanup(m);
		return 1;
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
			break;

		case fits:
			if( to_fits(m, fits_output) ) {
				perror("Error while saving image to FITS file");
				return 1;
			}
			break;

		case performance:
			printf("Ran %d iterations in %.3f [s] (%.3f [ms] per iteration)\n", iterations, duration/1000000., duration/1000./iterations);
			break;

		default:
			fprintf(stderr, "Output not currently supported: %d\n", output);
	}

	profit_cleanup(m);
	return 0;

}
