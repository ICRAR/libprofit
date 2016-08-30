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
#include <string>
#include <sstream>

#include "profit/coresersic.h"
#include "profit/ferrer.h"
#include "profit/king.h"
#include "profit/moffat.h"
#include "profit/profit.h"
#include "profit/psf.h"
#include "profit/sersic.h"
#include "profit/sky.h"

using namespace profit;
using namespace std;

class invalid_cmdline : exception {
public:
	invalid_cmdline(const string& what) : m_what(what) {}
	invalid_cmdline(const invalid_cmdline &e) : m_what(e.m_what) {}
	~invalid_cmdline() throw() {}
	const char *what() const throw() { return m_what.c_str(); }

private:
	std::string m_what;
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

#define _READ_DOUBLE_OR_FAIL(key, val, name, dst) \
	do { \
		char *endptr; \
		double tmp; \
		if ( key == name ) { \
			tmp = strtod(val.c_str(), &endptr); \
			if( tmp == 0 && endptr == val ) { \
				ostringstream os; \
				os << "Invalid double value for " << key << ": " << val; \
				throw invalid_cmdline(os.str()); \
			} \
			dst = tmp;\
			return true; \
		} \
	} while(0);

#define _READ_FROM_LONGINT_OR_FAIL(key, val, name, dst, TYPE) \
	do { \
		char *endptr; \
		long int tmp; \
		if ( key == name ) { \
			tmp = strtol(val.c_str(), &endptr, 10); \
			if( tmp == 0 && endptr == val ) { \
				ostringstream os; \
				os << "Invalid integer value for " << key << ": " << val; \
				throw invalid_cmdline(os.str()); \
			} \
			dst = (TYPE)tmp;\
			return true; \
		} \
	} while(0);

#define _READ_BOOL_OR_FAIL(key, val, name, dst) _READ_FROM_LONGINT_OR_FAIL(key, val, name, dst, bool)
#define _READ_UINT_OR_FAIL(key, val, name, dst) _READ_FROM_LONGINT_OR_FAIL(key, val, name, dst, unsigned int)


bool _keyval_to_sersic(Profile *p, const string &key, const string &val) {
	SersicProfile *s = static_cast<SersicProfile *>(p);
	_READ_DOUBLE_OR_FAIL(key, val, "xcen",  s->xcen);
	_READ_DOUBLE_OR_FAIL(key, val, "ycen",  s->ycen);
	_READ_DOUBLE_OR_FAIL(key, val, "mag",   s->mag);
	_READ_DOUBLE_OR_FAIL(key, val, "re",    s->re);
	_READ_DOUBLE_OR_FAIL(key, val, "nser",  s->nser);
	_READ_DOUBLE_OR_FAIL(key, val, "ang",   s->ang);
	_READ_DOUBLE_OR_FAIL(key, val, "axrat", s->axrat);
	_READ_DOUBLE_OR_FAIL(key, val, "box",   s->box);

	_READ_BOOL_OR_FAIL(  key, val, "rough",          s->rough);
	_READ_DOUBLE_OR_FAIL(key, val, "acc",            s->acc);
	_READ_DOUBLE_OR_FAIL(key, val, "rscale_switch",  s->rscale_switch);
	_READ_UINT_OR_FAIL(  key, val, "resolution",     s->resolution);
	_READ_UINT_OR_FAIL(  key, val, "max_recursions", s->max_recursions);
	_READ_BOOL_OR_FAIL(  key, val, "adjust",         s->adjust);

	_READ_DOUBLE_OR_FAIL(key, val, "rscale_max",     s->rscale_max);
	_READ_BOOL_OR_FAIL(  key, val, "rescale_flux",   s->rescale_flux);

	_READ_BOOL_OR_FAIL(key, val, "convolve", p->convolve);
	return false;
}

bool _keyval_to_moffat(Profile *p, const string &key, const string &val) {
	MoffatProfile *m = static_cast<MoffatProfile *>(p);
	_READ_DOUBLE_OR_FAIL(key, val, "xcen",  m->xcen);
	_READ_DOUBLE_OR_FAIL(key, val, "ycen",  m->ycen);
	_READ_DOUBLE_OR_FAIL(key, val, "mag",   m->mag);
	_READ_DOUBLE_OR_FAIL(key, val, "fwhm",  m->fwhm);
	_READ_DOUBLE_OR_FAIL(key, val, "con",   m->con);
	_READ_DOUBLE_OR_FAIL(key, val, "ang",   m->ang);
	_READ_DOUBLE_OR_FAIL(key, val, "axrat", m->axrat);
	_READ_DOUBLE_OR_FAIL(key, val, "box",   m->box);

	_READ_BOOL_OR_FAIL(  key, val, "rough",          m->rough);
	_READ_DOUBLE_OR_FAIL(key, val, "acc",            m->acc);
	_READ_DOUBLE_OR_FAIL(key, val, "rscale_switch",  m->rscale_switch);
	_READ_UINT_OR_FAIL(  key, val, "resolution",     m->resolution);
	_READ_UINT_OR_FAIL(  key, val, "max_recursions", m->max_recursions);
	_READ_BOOL_OR_FAIL(  key, val, "adjust",         m->adjust);

	_READ_DOUBLE_OR_FAIL(key, val, "rscale_max",     m->rscale_max);

	return false;
}

bool _keyval_to_ferrer(Profile *p, const string &key, const string &val) {
	FerrerProfile *f = static_cast<FerrerProfile *>(p);
	_READ_DOUBLE_OR_FAIL(key, val, "xcen",  f->xcen);
	_READ_DOUBLE_OR_FAIL(key, val, "ycen",  f->ycen);
	_READ_DOUBLE_OR_FAIL(key, val, "mag",   f->mag);
	_READ_DOUBLE_OR_FAIL(key, val, "rout",  f->rout);
	_READ_DOUBLE_OR_FAIL(key, val, "a",     f->a);
	_READ_DOUBLE_OR_FAIL(key, val, "b",     f->b);
	_READ_DOUBLE_OR_FAIL(key, val, "ang",   f->ang);
	_READ_DOUBLE_OR_FAIL(key, val, "axrat", f->axrat);
	_READ_DOUBLE_OR_FAIL(key, val, "box",   f->box);

	_READ_BOOL_OR_FAIL(  key, val, "rough",          f->rough);
	_READ_DOUBLE_OR_FAIL(key, val, "acc",            f->acc);
	_READ_DOUBLE_OR_FAIL(key, val, "rscale_switch",  f->rscale_switch);
	_READ_UINT_OR_FAIL(  key, val, "resolution",     f->resolution);
	_READ_UINT_OR_FAIL(  key, val, "max_recursions", f->max_recursions);
	_READ_BOOL_OR_FAIL(  key, val, "adjust",         f->adjust);

	_READ_DOUBLE_OR_FAIL(key, val, "rscale_max",     f->rscale_max);

	return false;
}

bool _keyval_to_coresersic(Profile *p, const string &key, const string &val) {
	CoreSersicProfile *cs = static_cast<CoreSersicProfile *>(p);
	_READ_DOUBLE_OR_FAIL(key, val, "xcen",  cs->xcen);
	_READ_DOUBLE_OR_FAIL(key, val, "ycen",  cs->ycen);
	_READ_DOUBLE_OR_FAIL(key, val, "mag",   cs->mag);
	_READ_DOUBLE_OR_FAIL(key, val, "re",    cs->re);
	_READ_DOUBLE_OR_FAIL(key, val, "nser",  cs->nser);
	_READ_DOUBLE_OR_FAIL(key, val, "rb",    cs->rb);
	_READ_DOUBLE_OR_FAIL(key, val, "a",     cs->a);
	_READ_DOUBLE_OR_FAIL(key, val, "b",     cs->b);
	_READ_DOUBLE_OR_FAIL(key, val, "ang",   cs->ang);
	_READ_DOUBLE_OR_FAIL(key, val, "axrat", cs->axrat);
	_READ_DOUBLE_OR_FAIL(key, val, "box",   cs->box);

	_READ_BOOL_OR_FAIL(  key, val, "rough",          cs->rough);
	_READ_DOUBLE_OR_FAIL(key, val, "acc",            cs->acc);
	_READ_DOUBLE_OR_FAIL(key, val, "rscale_switch",  cs->rscale_switch);
	_READ_UINT_OR_FAIL(  key, val, "resolution",     cs->resolution);
	_READ_UINT_OR_FAIL(  key, val, "max_recursions", cs->max_recursions);
	_READ_BOOL_OR_FAIL(  key, val, "adjust",         cs->adjust);

	_READ_DOUBLE_OR_FAIL(key, val, "rscale_max",     cs->rscale_max);

	_READ_BOOL_OR_FAIL(key, val, "convolve", p->convolve);
	return false;
}

bool _keyval_to_king(Profile *p, const string &key, const string &val) {
	KingProfile *cs = static_cast<KingProfile *>(p);
	_READ_DOUBLE_OR_FAIL(key, val, "xcen",  cs->xcen);
	_READ_DOUBLE_OR_FAIL(key, val, "ycen",  cs->ycen);
	_READ_DOUBLE_OR_FAIL(key, val, "mag",   cs->mag);
	_READ_DOUBLE_OR_FAIL(key, val, "rt",    cs->rt);
	_READ_DOUBLE_OR_FAIL(key, val, "rc",    cs->rc);
	_READ_DOUBLE_OR_FAIL(key, val, "a",     cs->a);
	_READ_DOUBLE_OR_FAIL(key, val, "ang",   cs->ang);
	_READ_DOUBLE_OR_FAIL(key, val, "axrat", cs->axrat);
	_READ_DOUBLE_OR_FAIL(key, val, "box",   cs->box);

	_READ_BOOL_OR_FAIL(  key, val, "rough",          cs->rough);
	_READ_DOUBLE_OR_FAIL(key, val, "acc",            cs->acc);
	_READ_DOUBLE_OR_FAIL(key, val, "rscale_switch",  cs->rscale_switch);
	_READ_UINT_OR_FAIL(  key, val, "resolution",     cs->resolution);
	_READ_UINT_OR_FAIL(  key, val, "max_recursions", cs->max_recursions);
	_READ_BOOL_OR_FAIL(  key, val, "adjust",         cs->adjust);

	_READ_DOUBLE_OR_FAIL(key, val, "rscale_max",     cs->rscale_max);

	_READ_BOOL_OR_FAIL(key, val, "convolve", p->convolve);
	return false;
}

bool _keyval_to_sky(Profile *p, const string &key, const string &val) {
	SkyProfile *s = static_cast<SkyProfile *>(p);
	_READ_DOUBLE_OR_FAIL(key, val, "bg",  s->bg);
	_READ_BOOL_OR_FAIL(key, val, "convolve", p->convolve);
	return false;
}

bool _keyval_to_psf(Profile *p, const string &key, const string &val) {
	PsfProfile *s = static_cast<PsfProfile *>(p);
	_READ_DOUBLE_OR_FAIL(key, val, "xcen",  s->xcen);
	_READ_DOUBLE_OR_FAIL(key, val, "ycen",  s->ycen);
	_READ_DOUBLE_OR_FAIL(key, val, "mag",   s->mag);
	return false;
}

void desc_to_profile(
	Model &model,
	string description,
	const char* name,
	bool (keyval_to_param)(Profile *, const string& name, const string &value)
) {

	string tok;
	Profile *p;
	bool assigned;

	p = model.add_profile(name);
	if( p == NULL ) {
		ostringstream os;
		os << "No profile found for profile name: " << name;
		throw invalid_cmdline(os.str());
	}
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

		assigned = keyval_to_param(p, name_and_value[0], name_and_value[1]);
		if( !assigned ) {
			cerr << "Ignoring unknown " << name << " profile parameter: " << name_and_value[0] << endl;
		}

	}
}

void parse_profile(Model &model, const string &description) {

	/* The description might be only a name */
	string subdesc;
	string::size_type colon = description.find(':');
	if( colon != string::npos ) {
		subdesc = description.substr(colon + 1);
	}

	if( !description.compare(0, 6, "sersic") ) {
		desc_to_profile(model, subdesc, "sersic", &_keyval_to_sersic);
	}
	else if( !description.compare(0, 6, "moffat") ) {
		desc_to_profile(model, subdesc, "moffat", &_keyval_to_moffat);
	}
	else if( !description.compare(0, 6, "ferrer") ) {
		desc_to_profile(model, subdesc, "ferrer", &_keyval_to_ferrer);
	}
	else if( !description.compare(0, 10, "coresersic") ) {
		desc_to_profile(model, subdesc, "coresersic", &_keyval_to_coresersic);
	}
	else if( !description.compare(0, 4, "king") ) {
		desc_to_profile(model, subdesc, "king", &_keyval_to_king);
	}
	else if( !description.compare(0, 3, "sky") ) {
		desc_to_profile(model, subdesc, "sky", &_keyval_to_sky);
	}
	else if( !description.compare(0, 3, "psf") ) {
		desc_to_profile(model, subdesc, "psf", &_keyval_to_psf);
	}
	else {
		ostringstream os;
		os << "Unknown profile name in profile description: " << description;
		throw invalid_cmdline(os.str());
	}
}

double *parse_psf(string optarg,
                  unsigned int &psf_width, unsigned int &psf_height,
                  double &psf_scale_x, double &psf_scale_y) {

	char *endptr;
	const char *tok;
	unsigned int size, i = 0;
	double *psf;
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
	tok = (*it++).c_str();
	psf_width = (unsigned int)strtoul(tok, &endptr, 10);
	if( tok == endptr ) {
		ostringstream os;
		os << "Invalid value for psf's width: " << tok;
		throw invalid_cmdline(os.str());
	}

	tok = (*it++).c_str();
	psf_height = (unsigned int)strtoul(tok, &endptr, 10);
	if( tok == endptr ) {
		ostringstream os;
		os << "Invalid value for psf's height: " << tok;
		throw invalid_cmdline(os.str());
	}

	if( read_scales ) {
		tok = (*it++).c_str();
		psf_scale_x = strtod(tok, &endptr);
		if( tok == endptr ) {
			ostringstream os;
			os << "Invalid value for psf's scale_x: " << tok;
			throw invalid_cmdline(os.str());
		}

		tok = (*it++).c_str();
		psf_scale_y = strtod(tok, &endptr);
		if( tok == endptr ) {
			ostringstream os;
			os << "Invalid value for psf's scale_y: " << tok;
			throw invalid_cmdline(os.str());
		}
	}

	size = psf_width * psf_height;
	psf = new double[size];

	vector<string> values;
	tokenize(*it, values, ",");
	i = 0;
	for(auto value: values) {
		const char *v = value.c_str();
		psf[i] = strtod(v, &endptr);
		if( psf[i] == 0 && v == endptr ) {
			delete [] psf;
			ostringstream os;
			os << "Invalid floating-point value for psf: " << value;
			throw invalid_cmdline(os.str());
		}
		i++;
	}

	if( i != size ) {
		delete [] psf;
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

double *read_image_from_fits_file(const string &filename, unsigned int &width, unsigned int &height, double &scale_x, double &scale_y) {

	FILE *f;
	unsigned int i, pos, padding;
	char hdr[80];

	width = height = 0;
	scale_x = scale_y = 1;

	f = fopen(filename.c_str(), "rb");
	if( !f ) {
		perror("Couldn't open file for reading");
		return NULL;
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
	double *out = new double[size];
	size_t nitems = fread(out, sizeof(double), size, f);
	fclose(f);
	if( nitems != size ) {
		perror("Error while reading file");
		delete [] out;
		return NULL;
	}

	/* data has to be big-endian */
	if( is_little_endian() ) {
		double *it = out;
		for(i=0; i!=size; i++) {
			*it = swap_bytes(*it);
			it++;
		}
	}

	return out;
}

int to_fits(Model &m, vector<double> image, string fname) {

	FILE *f;
	unsigned int i, j, pos, padding;
	char hdr[80];

	/* Append .fits if not in the name yet */
	size_t fname_size = fname.size();
	if( fname_size <= 5 || fname.rfind(".fits", fname_size - 6) != string::npos ) {
		stringstream ss;
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
				if( !m.psf ) {
					usage(stderr, argv);
					return 1;
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
				if( magzero == 0 && endptr == optarg ) {
					ostringstream os;
					os << "Invalid magzero value: " << optarg;
					throw invalid_cmdline(os.str());
				}
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
	if( !m.profiles.size() ) {
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
