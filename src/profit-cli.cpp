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
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <sstream>

#include "profit/profit.h"


using namespace std;
using namespace profit;

class invalid_cmdline : exception {
public:
	invalid_cmdline(const string& what) : m_what(what) {}
	invalid_cmdline(const invalid_cmdline &e) : m_what(e.m_what) {}
	~invalid_cmdline() throw() {}
	const char *what() const throw() { return m_what.c_str(); }

private:
	string m_what;
};

/**
 * Breaks down a string into substrings delimited by delims
 *
 * Taken from:
 *  http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
 */
void tokenize(const string &s, vector<string> &tokens, const string &delims) {

	string::size_type lastPos = s.find_first_not_of(delims, 0);
	string::size_type pos     = s.find_first_of(delims, lastPos);

	while (string::npos != pos || string::npos != lastPos) {
		tokens.push_back(s.substr(lastPos, pos - lastPos));
		lastPos = s.find_first_not_of(delims, pos);
		pos = s.find_first_of(delims, lastPos);
	}

}

template <typename T, typename F>
void read_from_string(Profile &p, const string& key, string &val, const string &name, F reader) {

	if( val.empty() || key != name ) {
		return;
	}

	try {
		T tgt = reader(val);
		p.parameter(name, tgt);
		val.erase();
	} catch(invalid_argument &e) {
		ostringstream os;
		os << "Invalid double value '" << val << "' for " << key;
		throw invalid_cmdline(os.str());
	}
}

static
void read_dble(Profile &p, const string &key, string &val, const string &name) {
	read_from_string<double>(p, key, val, name, [](const string &v){return stod(v);});
}

static
void read_uint(Profile &p, const string &key, string &val, const string &name) {
	read_from_string<unsigned int>(p, key, val, name, [](const string &v){return stoul(v, nullptr, 10);});
}

static
void read_bool(Profile &p, const string &key, string &val, const string &name) {
	read_from_string<bool>(p, key, val, name, [](const string &v){return stoul(v, nullptr, 10);});
}

static
void _keyval_to_radial(Profile &p, const string &key, string &val) {
	read_dble(p, key, val, "xcen");
	read_dble(p, key, val, "ycen");
	read_dble(p, key, val, "mag");
	read_dble(p, key, val, "ang");
	read_dble(p, key, val, "axrat");
	read_dble(p, key, val, "box");

	read_bool(p, key, val, "rough");
	read_dble(p, key, val, "acc");
	read_dble(p, key, val, "rscale_switch");
	read_uint(p, key, val, "resolution");
	read_uint(p, key, val, "max_recursions");
	read_bool(p, key, val, "adjust");
	read_dble(p, key, val, "rscale_max");
}

static
void keyval_to_sersic(Profile &p, const string &key, string &val) {
	_keyval_to_radial(p, key, val);
	read_dble(p, key, val, "re");
	read_dble(p, key, val, "nser");
	read_bool(p, key, val, "rescale_flux");
}

static
void keyval_to_moffat(Profile &p, const string &key, string &val) {
	_keyval_to_radial(p, key, val);
	read_dble(p, key, val, "fwhm");
	read_dble(p, key, val, "con");
}

static
void keyval_to_ferrer(Profile &p, const string &key, string &val) {
	_keyval_to_radial(p, key, val);
	read_dble(p, key, val, "rout");
	read_dble(p, key, val, "a");
	read_dble(p, key, val, "b");
}

static
void keyval_to_coresersic(Profile &p, const string &key, string &val) {
	_keyval_to_radial(p, key, val);
	read_dble(p, key, val, "re");
	read_dble(p, key, val, "nser");
	read_dble(p, key, val, "rb");
	read_dble(p, key, val, "a");
	read_dble(p, key, val, "b");
}

static
void keyval_to_brokenexp(Profile &p, const string &key, string &val) {
	_keyval_to_radial(p, key, val);
	read_dble(p, key, val, "h1");
	read_dble(p, key, val, "h2");
	read_dble(p, key, val, "rb");
	read_dble(p, key, val, "a");
}

static
void keyval_to_king(Profile &p, const string &key, string &val) {
	_keyval_to_radial(p, key, val);
	read_dble(p, key, val, "rt");
	read_dble(p, key, val, "rc");
	read_dble(p, key, val, "a");
}

static
void keyval_to_sky(Profile &p, const string &key, string &val) {
	read_dble(p, key, val, "bg");
}

static
void keyval_to_psf(Profile &p, const string &key, string &val) {
	read_dble(p, key, val, "xcen");
	read_dble(p, key, val, "ycen");
	read_dble(p, key, val, "mag");
}

typedef void (*keyval_to_param_t)(Profile &, const string& name, string &value);
static map<string, keyval_to_param_t> reader_functions = {
	{"sersic",     &keyval_to_sersic},
	{"moffat",     &keyval_to_moffat},
	{"ferrer",     &keyval_to_ferrer},
	{"ferrers",    &keyval_to_ferrer},
	{"king",       &keyval_to_king},
	{"coresersic", &keyval_to_coresersic},
	{"brokenexp",  &keyval_to_brokenexp},
	{"sky",        &keyval_to_sky},
	{"psf",        &keyval_to_psf}
};

void desc_to_profile(
	Model &model,
	const string &name,
	string description
) {

	string tok;

	Profile &p = model.add_profile(name);
	keyval_to_param_t keyval_to_param = reader_functions[name];

	if( description.size() == 0 ) {
		return;
	}

	vector<string> tokens;
	tokenize(description, tokens, ":");
	for(auto token: tokens) {

		vector<string> name_and_value;
		tokenize(token, name_and_value, "=");
		if( name_and_value.size() != 2 ) {
			ostringstream os;
			os <<  "Parameter " << token << " of profile " << name << " doesn't obey the form name=value";
			throw invalid_cmdline(os.str());
		}

		const string &key = name_and_value[0];
		string &val = name_and_value[1];
		read_bool(p, key, val, "convolve");
		keyval_to_param(p, key, val);
		if( !name_and_value[1].empty() ) {
			cerr << "Ignoring unknown " << name << " profile parameter: " << name_and_value[0] << endl;
		}

	}
}

void parse_profile(Model &model, const string &description) {

	/* The description might be only a name */
	string name;
	string subdesc;
	string::size_type colon = description.find(':');
	if( colon != string::npos ) {
		name = description.substr(0, colon);
		subdesc = description.substr(colon + 1);
	}
	else {
		name = description;
	}

	desc_to_profile(model, name, subdesc);
}

vector<double> parse_psf(string optarg,
                         unsigned int &psf_width, unsigned int &psf_height,
                         double &psf_scale_x, double &psf_scale_y) {

	unsigned int size, i = 0;
	bool read_scales = false;

	/* format is w:h:[optional scale_x:scale_y:]:val1,val2... */
	vector<string> tokens;
	tokenize(optarg, tokens, ":");
	vector<string>::size_type ntokens = tokens.size();
	if( ntokens < 1 ) {
		throw invalid_cmdline("Missing psf's width");
	}
	else if( ntokens < 2 ) {
		throw invalid_cmdline("Missing psf's height");
	}
	else if( ntokens != 3 && ntokens != 5 ) {
		throw invalid_cmdline("Invalid psf format, see -h for help");
	}
	else if( ntokens == 5 ) {
		read_scales = true;
	}

	vector<string>::iterator it = tokens.begin();
	psf_width = stoul(*it++);
	psf_height = stoul(*it++);
	if( read_scales ) {
		psf_scale_x = stod(*it++);
		psf_scale_y = stod(*it++);
	}

	size = psf_width * psf_height;
	vector<double> psf(size);

	vector<string> values;
	tokenize(*it, values, ",");
	i = 0;
	for(auto value: values) {
		psf[i++] = stod(value);
	}

	if( i != size ) {
		ostringstream os;
		os << "Not enough values provided for PSF. Provided: " << i << ", expected: " << size;
		throw invalid_cmdline(os.str());
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
	fprintf(file,"  -x        Image width. Defaults to 100\n");
	fprintf(file,"  -y        Image height. Defaults to 100\n");
	fprintf(file,"  -w        Width in pixels. Defaults to 100\n");
	fprintf(file,"  -H        Height in pixels. Defaults to 100\n");
	fprintf(file,"  -m        Zero magnitude. Defaults to 0.\n");
	fprintf(file,"  -P        PSF function (specified as w:h:val1,val2..., or as a FITS filename)\n");
	fprintf(file,"  -h,-?     Show this help and exit\n");
	fprintf(file,"  -v        Show the program version and exit\n\n");
	fprintf(file,"Profiles should be specified as follows:\n\n");
	fprintf(file,"-p name:param1=val1:param2=val2:...\n\n");
	fprintf(file,"The following profiles (and parameters) are currently accepted:\n\n");
	fprintf(file," * psf: xcen, ycen, mag\n");
	fprintf(file," * sky: bg\n");
	fprintf(file," * sersic: re, nser, rescale_flux\n");
	fprintf(file," * moffat: fwhm, con\n");
	fprintf(file," * ferrer: a, b, rout\n");
	fprintf(file," * coresersic: re, nser, rb, a, b\n");
	fprintf(file," * brokenexp: h1, h2, rb, a\n");
	fprintf(file," * king: rc, rt, a\n");
	fprintf(file,"\
 * sersic, moffat, ferrer, coresersic, king: xcen, ycen, mag, box, ang, axrat,\n\
                           rough, rscale_switch, max_recursions,\n\
                           resolution, acc, rscale_max, adjust\n\n");
	fprintf(file,"For more information visit https://libprofit.readthedocs.io.\n\n");
}

static inline
bool is_little_endian() {
	volatile uint32_t i=0x01234567;
	return (*((uint8_t*)(&i))) == 0x67;
}

static inline
double swap_bytes(double v) {
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

#define FITS_BLOCK_SIZE (36*80)

vector<double> read_image_from_fits_file(const string &filename, unsigned int &width, unsigned int &height, double &scale_x, double &scale_y) {

	FILE *f;
	unsigned int pos, padding;
	char hdr[80];

	width = height = 0;
	scale_x = scale_y = 1;

	f = fopen(filename.c_str(), "rb");
	if( !f ) {
		ostringstream ss;
		ss << "Couldn't open '" << filename << "' for reading: " << strerror(errno);
		throw invalid_cmdline(ss.str());
	}

	/*
	 * Standard headers, we're assuming they say 'T" for SIMPLE, -64 for BITPIX
	 * and 2 for NAXIS.
	 */
	while( fread(hdr, 1, 80, f) ) {

		if( !strncmp("NAXIS1", hdr, 6) ) {
			sscanf(hdr, "NAXIS1 = %u", &width);
		}
		else if( !strncmp("NAXIS2", hdr, 6) ) {
			sscanf(hdr, "NAXIS2 = %u", &height);
		}
		else if( !strncmp("CDELT1", hdr, 6) ) {
			sscanf(hdr, "CDELT1 = %lf", &scale_x);
		}
		else if( !strncmp("CDELT2", hdr, 6) ) {
			sscanf(hdr, "CDELT2 = %lf", &scale_y);
		}
		else if( !strncmp("END", hdr, 3) ) {
			break;
		}
	}

	pos = (unsigned int)ftell(f);
	padding = FITS_BLOCK_SIZE - (pos % FITS_BLOCK_SIZE);
	fseek(f, padding, SEEK_CUR);

	unsigned int size = width * height;
	vector<double> psf(size);
	size_t nitems = fread(psf.data(), sizeof(double), size, f);
	fclose(f);
	if( nitems != size ) {
		throw invalid_cmdline("Error while reading file: less data found than expected");
	}

	/* data has to be big-endian */
	if( is_little_endian() ) {
		transform(psf.begin(), psf.end(), psf.begin(), swap_bytes);
	}

	return psf;
}

int to_fits(Model &m, vector<double> image, string fname) {

	FILE *f;
	unsigned int i, pos, padding;
	char hdr[80];

	/* Append .fits if not in the name yet */
	size_t fname_size = fname.size();
	if( fname_size <= 5 || fname.rfind(".fits", fname_size - 6) != string::npos ) {
		ostringstream ss;
		ss << fname << ".fits";
		fname = ss.str();
	}

	f = fopen(fname.c_str(), "w+");
	if( !f ) {
		return 1;
	}

	/*
	 * Standard headers
	 *
	 * The first five headers are required, and must be in "fixed format",
	 * meaning that their values must be right-indented on column 30, sigh...
	 */
	fprintf(f, "%-80s", "SIMPLE  =                    T / File conforms to FITS standard");
	fprintf(f, "%-80s", "BITPIX  =                  -64 / Bits per pixel");
	fprintf(f, "%-80s", "NAXIS   =                    2 / Number of axes");
	sprintf(hdr, "NAXIS1  =           %10.0u / Width", m.width);
	fprintf(f, "%-80s", hdr);
	sprintf(hdr, "NAXIS2  =           %10.0u / Height", m.height);
	fprintf(f, "%-80s", hdr);
	fprintf(f, "%-80s", "CRPIX1  = 1");
	sprintf(hdr, "CRVAL1  = %f", 0.5*m.scale_x);
	fprintf(f, "%-80s", hdr);
	sprintf(hdr, "CDELT1  = %f", m.scale_x);
	fprintf(f, "%-80s", hdr);
	fprintf(f, "%-80s", "CTYPE1  = ' '");
	fprintf(f, "%-80s", "CUNIT1  = ' '");
	fprintf(f, "%-80s", "CRPIX2  = 1");
	sprintf(hdr, "CRVAL2  = %f", 0.5*m.scale_y);
	fprintf(f, "%-80s", hdr);
	sprintf(hdr, "CDELT2  = %f", m.scale_y);
	fprintf(f, "%-80s", hdr);
	fprintf(f, "%-80s", "CTYPE2  = ' '");
	fprintf(f, "%-80s", "CUNIT2  = ' '");
	fprintf(f, "%-80s", "END");

	pos = (unsigned int)ftell(f);
	padding = FITS_BLOCK_SIZE - (pos % FITS_BLOCK_SIZE);
	for(i=0; i<padding; i++) {
		fprintf(f, " ");
	}

	/* data has to be big-endian */
	size_t image_size = m.width * m.height;
	if( is_little_endian() ) {
		vector<double> big_endian_image(image_size);
		transform(image.begin(), image.end(), big_endian_image.begin(), swap_bytes);
		fwrite(big_endian_image.data(), sizeof(double), image_size, f);
	}
	else {
		fwrite(image.data(), sizeof(double), image_size, f);
	}

	/* Pad with zeroes until we complete the current 36*80 block */
	padding = FITS_BLOCK_SIZE - (((unsigned int)sizeof(double) * image_size) % FITS_BLOCK_SIZE);
	string zeros(padding, 0);
	fwrite(zeros.c_str(), 1, padding, f);
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

int parse_and_run(int argc, char *argv[]) {
	int opt;
	unsigned int width = 100, height = 100, iterations = 1;
	long duration;
	double magzero = 0, scale_x = 1, scale_y = 1;
	unsigned int i, j;
	char *endptr = NULL;
	string fits_output;
	output_t output = none;
	Model m;
	struct stat stat_buf;

	while( (opt = getopt(argc, argv, "h?vP:p:w:H:x:y:X:Y:m:tbf:i:")) != -1 ) {
		switch(opt) {

			case 'h':
			case '?':
				usage(stdout, argv);
				return 0;

			case 'v':
				printf("libprofit version %s\n", PROFIT_VERSION);
				return 0;

			case 'p':
				parse_profile(m, optarg);
				break;

			case 'P':
				if( !stat(optarg, &stat_buf) ) {
					m.psf = read_image_from_fits_file(optarg, m.psf_width, m.psf_height, m.psf_scale_x, m.psf_scale_y);
				}
				else {
					m.psf = parse_psf(optarg, m.psf_width, m.psf_height, m.psf_scale_x, m.psf_scale_y);
				}
				break;

			case 'w':
				width = (unsigned int)atoi(optarg);
				break;

			case 'H':
				height = (unsigned int)atoi(optarg);
				break;

			case 'x':
				scale_x = atof(optarg);
				break;

			case 'y':
				scale_y = atof(optarg);
				break;

			case 'm':
				magzero = strtod(optarg, &endptr);
				break;

			case 't':
				if( output != none ) {
					throw invalid_cmdline("-t and -b cannot be used together");
				}
				output = text;
				break;

			case 'b':
				if( output != none ) {
					throw invalid_cmdline("-b and -t cannot be used together");
				}
				output = binary;
				break;

			case 'f':
				fits_output = optarg;
				output = fits;
				break;

			case 'i':
				iterations = (unsigned int)atoi(optarg);
				output = performance;
				break;

			default:
				usage(stderr, argv);
				return 1;

		}
	}

	/* No profiles given */
	if( !m.has_profiles() ) {
		usage(stderr, argv);
		return 1;
	}

	/* We default to text output */
	if( output == none ) {
		output = text;
	}

	m.width   = width;
	m.height  = height;
	m.scale_x = scale_x;
	m.scale_y = scale_y;
	m.magzero = magzero;

	/* This means that we evaluated the model once, but who cares */
	struct timeval start, end;
	gettimeofday(&start, NULL);
	vector<double> image;
	for(i=0; i!=iterations; i++) {
		image = m.evaluate();
	}
	gettimeofday(&end, NULL);
	duration = (end.tv_sec - start.tv_sec)*1000000 + (end.tv_usec - start.tv_usec);

	switch(output) {

		case binary:
			fwrite(image.data(), sizeof(double), m.width * m.height, stdout);
			break;

		case text:
			for(j=0; j!=m.height; j++) {
				for(i=0; i!=m.width; i++) {
					printf("%g ", image[j*m.width + i]);
				}
				printf("\n");
			}
			break;

		case fits:
			if( to_fits(m, image, fits_output) ) {
				perror("Error while saving image to FITS file");
				return 1;
			}
			break;

		case performance:
			printf("Ran %d iterations in %.3f [s] (%.3f [ms] per iteration)\n", iterations, (double)duration/1000000., (double)duration/1000./iterations);
			break;

		default:
			cerr << "Output not currently supported: " << output << endl;
	}

	return 0;

}

int main(int argc, char *argv[]) {
	try {
		return parse_and_run(argc, argv);
	} catch (invalid_cmdline &e) {
		cerr << "Error on command line: " << e.what() << endl;
		return 1;
	} catch (invalid_parameter &e) {
		cerr << "Error while calculating model: " << e.what() << endl;
		return 1;
	}
}
