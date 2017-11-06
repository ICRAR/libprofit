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
#include <sys/types.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
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

	shared_ptr<Profile> p = model.add_profile(name);
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
		read_bool(*p, key, val, "convolve");
		keyval_to_param(*p, key, val);
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
	fprintf(file,"  -s        Show runtime stats\n");
	fprintf(file,"  -T <conv> Use this type of convolver (see below)\n");
#ifdef PROFIT_OPENCL
	fprintf(file,"  -C <p,d>  Use OpenCL with platform p, device d, and double support (0|1)\n");
	fprintf(file,"  -c        Display OpenCL information about devices and platforms\n");
#endif /* PROFIT_OPENCL */
#ifdef PROFIT_OPENMP
	fprintf(file,"  -n <n>    Use n OpenMP threads to calculate profiles\n");
#endif /* PROFIT_OPENMP */
#ifdef PROFIT_FFTW
	fprintf(file,"  -F <n>    FFTW plans created with n effort (more takes longer)\n");
	fprintf(file,"  -r        Reuse FFT-transformed PSF across evaluations (if -T fft)\n");
#endif /* PROFIT_FFTW */
	fprintf(file,"  -x        Image width. Defaults to 100\n");
	fprintf(file,"  -y        Image height. Defaults to 100\n");
	fprintf(file,"  -w        Width in pixels. Defaults to 100\n");
	fprintf(file,"  -H        Height in pixels. Defaults to 100\n");
	fprintf(file,"  -m        Zero magnitude. Defaults to 0.\n");
	fprintf(file,"  -P        PSF function (specified as w:h:val1,val2..., or as a FITS filename)\n");
	fprintf(file,"  -h,-?     Show this help and exit\n");
	fprintf(file,"  -V        Show the program version and exit\n\n");
	fprintf(file,"The following convolver types are supported:\n\n");
	fprintf(file," * brute: A brute-force convolver\n");
	fprintf(file," * brute-old: An older, slower brute-force convolver (used only for comparisons)\n");
#ifdef PROFIT_OPENCL
	fprintf(file," * opencl: An OpenCL-based brute-force convolver\n");
	fprintf(file," * opencl-local: An OpenCL-based local-memory caching brute-force convolver\n");
#endif // PROFIT_OPENCL
#ifdef PROFIT_FFTW
	fprintf(file," * fft: An FFT-based convolver\n");
#endif // PROFIT_FFTW
	fprintf(file,"\nProfiles should be specified as follows:\n\n");
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

static
void print_stats_line(const string &prefix, const string &stat_name, double val) {
	const auto static name_width = 50u;
	int nchars = prefix.size() + stat_name.size();
	int nspaces = name_width - nchars;
	string spaces(max(0, nspaces), ' ');
	cout << prefix << stat_name << spaces << " : " << setw(10) << setprecision(3) << setiosflags(ios::fixed) << val << " [ms]" << endl;
}

#ifdef PROFIT_OPENCL
static
void print_opencl_info() {

	const auto info = get_opencl_info();

	if( info.size() > 0 ) {
		cout << "OpenCL information" << endl;
		cout << "==================" << endl << endl;
		for(auto platform_info: info) {
			auto plat_id = get<0>(platform_info);
			auto plat_info = get<1>(platform_info);
			cout << "Platform [" << plat_id << "]" << endl;
			cout << "  Name           : " << plat_info.name << endl;
			cout << "  OpenCL version : " << plat_info.supported_opencl_version/100. << endl;
			for(auto device_info: plat_info.dev_info) {
				cout << "  Device [" << get<0>(device_info) << "]" << endl;
				cout << "    Name         : " << get<1>(device_info).name << endl;
				cout << "    Double       : " << (get<1>(device_info).double_support ? "Supported" : "Not supported") << endl;
			}
			cout << endl;
		}
	}
	else {
		cout << "No OpenCL installation found" << endl;
	}
}

static
void print_cl_stats(const string &prefix0, bool opencl_120, const OpenCL_times &stats) {

	auto prefix1 = prefix0 + "  ";

	ostringstream os;
	os << "OpenCL operations (" << stats.nwork_items << " work items)";
	print_stats_line(prefix0, os.str(), stats.total / 1e6 );
	print_stats_line(prefix1, "Kernel preparation", stats.kernel_prep / 1e6 );
	if( opencl_120 ) {
		print_stats_line(prefix1, "Fill submission", stats.filling_times.submit / 1e6 );
		print_stats_line(prefix1, "Fill execution", stats.filling_times.exec / 1e6 );
	}
	print_stats_line(prefix1, "Write submission", stats.writing_times.submit / 1e6 );
	print_stats_line(prefix1, "Write execution", stats.writing_times.exec / 1e6 );
	print_stats_line(prefix1, "Kernel submission", stats.kernel_times.submit / 1e6 );
	print_stats_line(prefix1, "Kernel execution", stats.kernel_times.exec / 1e6 );
	print_stats_line(prefix1, "Read submission", stats.reading_times.submit / 1e6 );
	print_stats_line(prefix1, "Read execution", stats.reading_times.exec / 1e6 );
}
#endif /* PROFIT_OPENCL */

static
void print_stats(const Model &m) {

#ifdef PROFIT_DEBUG
	for(const auto &profile_integrations: m.get_profile_integrations()) {
		int total = 0;
		if( get<1>(profile_integrations).size() > 0 ) {
			cout << "Integrations per recursion level for profile " << get<0>(profile_integrations) << endl;
			for(const auto level_integrations: get<1>(profile_integrations)) {
				auto integrations = get<1>(level_integrations);
				total += integrations;
				cout << " Level " << get<0>(level_integrations) << ": " << integrations << " integrations" << endl;
			}
			cout << " Total: " << total << " integrations" << endl;
		}
		else {
			cout << "Profile " << get<0>(profile_integrations) << " didn't run into any recursion" << endl;
		}
	}
#endif /* PROFIT_DEBUG */

	cout << endl;
	auto const &stats = m.get_stats();

	auto prefix0 = "";
	for(auto const &stat_pair: stats) {

		// Some profile might not have gathered stats
		auto profile_name = get<0>(stat_pair);
		auto stat = get<1>(stat_pair);
		ProfileStats *profile_stats = stat.get();
		if( !profile_stats ) {
			continue;
		}

		cout << "Stats for profile " << profile_name << endl;

#ifdef PROFIT_OPENCL
		auto prefix1 = "  ";
		RadialProfileStats *rprofile_stats = dynamic_cast<RadialProfileStats *>(profile_stats);
		if( rprofile_stats && m.opencl_env ) {
			bool opencl_120 = m.opencl_env->get_version() >= 120;
			print_cl_stats(prefix0, opencl_120, rprofile_stats->cl_times);
			print_stats_line(prefix0, "Pre-loop", rprofile_stats->subsampling.pre_subsampling / 1e6 );
			print_stats_line(prefix0, "Subsampling loop", rprofile_stats->subsampling.total / 1e6 );
			print_stats_line(prefix1, "New subsamples calculation", rprofile_stats->subsampling.new_subsampling / 1e6 );
			print_stats_line(prefix1, "Initial transform", rprofile_stats->subsampling.inital_transform / 1e6 );
			print_cl_stats(prefix1, opencl_120, rprofile_stats->subsampling.cl_times);
			print_stats_line(prefix1, "Final transform", rprofile_stats->subsampling.final_transform / 1e6 );
			print_stats_line(prefix0, "Final image", rprofile_stats->final_image / 1e6 );
		}
#endif /* PROFIT_OPENCL */

		print_stats_line(prefix0, "Total", profile_stats->total / 1e6 );
	}
}

static inline
bool is_little_endian() {
	const uint32_t i = 0x01234567;
	return *reinterpret_cast<const uint8_t *>(&i) == 0x67;
}

static inline
double swap_bytes(const double v) {
	double r;
	const char *vbytes = reinterpret_cast<const char *>(&v);
	char *rbytes = reinterpret_cast<char *>(&r);
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

int to_fits(Model &m, const Image image, string fname) {

	FILE *f;
	unsigned int i, pos, padding;
	char hdr[80];

	/* Append .fits if not in the name yet */
	size_t fname_size = fname.size();
	if( fname_size <= 5 || fname.compare(fname_size - 5, fname_size, ".fits") != 0 ) {
		fname = fname + ".fits";
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
	sprintf(hdr, "NAXIS1  =           %10.0u / Width", image.getWidth());
	fprintf(f, "%-80s", hdr);
	sprintf(hdr, "NAXIS2  =           %10.0u / Height", image.getHeight());
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
	size_t image_size = image.size();
	auto &imdata = image.getData();
	if( is_little_endian() ) {
		vector<double> big_endian_image(image_size);
		transform(imdata.begin(), imdata.end(), big_endian_image.begin(), swap_bytes);
		fwrite(big_endian_image.data(), sizeof(double), image_size, f);
	}
	else {
		fwrite(imdata.data(), sizeof(double), image_size, f);
	}

	/* Pad with zeroes until we complete the current 36*80 block */
	padding = FITS_BLOCK_SIZE - (((unsigned int)sizeof(double) * image_size) % FITS_BLOCK_SIZE);
	string zeros(padding, 0);
	fwrite(zeros.c_str(), 1, padding, f);
	fclose(f);

	return 0;
}

Image run(unsigned int iterations, Model &m) {

	using chrono::system_clock;

	/* This means that we evaluated the model once, but who cares */
	Image image;
	auto start = system_clock::now();
	for(unsigned i=0; i!=iterations; i++) {
		image = m.evaluate();
	}
	auto end = system_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(end-start).count();

	double dur_secs = (double)duration/1000;
	double dur_per_iter = (double)duration/iterations;
	cout << std::fixed << std::setprecision(3);
	cout << "Ran " << iterations << " iterations in ";
	cout << setprecision(3) << fixed << dur_secs << " [s] ";
	cout << "(" << setprecision(3) << fixed << dur_per_iter << " [ms] per iteration)";
	cout << endl;

	return image;
}

typedef enum _output_type {
	none = 0,
	binary = 1,
	text = 2,
	fits = 3,
} output_t;

int parse_and_run(int argc, char *argv[]) {

	using namespace std::chrono;
	using chrono::system_clock;

	int opt;
	unsigned int width = 100, height = 100, iterations = 1;
	double magzero = 0, scale_x = 1, scale_y = 1;
	unsigned int i, j;
	char *endptr = NULL;
	string fits_output;
	output_t output = none;
	Model m;
	string convolver_type = "brute";
	ConvolverCreationPreferences convolver_prefs;
	struct stat stat_buf;
	bool show_stats = false;

#ifdef PROFIT_OPENCL
	bool use_opencl = false, use_double = false;
	unsigned int clplat_idx = 0, cldev_idx = 0;
	vector<string> tokens;
#endif /* PROFIT_OPENCL */

	const char *options = "h?VsP:p:w:H:x:y:X:Y:m:tbf:i:T:"
#ifdef PROFIT_OPENCL
	                      "C:c"
#endif /* PROFIT_OPENCL */
#ifdef PROFIT_OPENMP
	                      "n:"
#endif /* PROFIT_OPENMP */
#ifdef PROFIT_FFTW
	                      "F:r"
#endif /* PROFIT_FFTW */
	;

	while( (opt = getopt(argc, argv, options)) != -1 ) {
		switch(opt) {

			case 'h':
			case '?':
				usage(stdout, argv);
				return 0;

			case 'V':
				cout << "libprofit version " << PROFIT_VERSION_MAJOR << "."
				     << PROFIT_VERSION_MINOR << "." << PROFIT_VERSION_PATCH << endl;
				cout << "OpenCL support: ";
#ifdef PROFIT_OPENCL
				cout << "Yes (up to " << PROFIT_OPENCL_MAJOR << "." << PROFIT_OPENCL_MINOR << ")";
#else
				cout << "No";
#endif
				cout << endl << "OpenMP support: ";
#ifdef PROFIT_OPENMP
				cout << "Yes";
#else
				cout << "No";
#endif
				cout << endl;
				cout << "FFTW support: ";
#ifdef PROFIT_FFTW
				cout << "Yes";
#else
				cout << "No";
#endif /* PROFIT_FFTW */
				cout << endl;
				return 0;

			case 's':
				show_stats = true;
				break;

			case 'T':
				convolver_type = optarg;
				break;

#ifdef PROFIT_FFTW
			case 'F':
				convolver_prefs.effort = FFTPlan::effort_t(std::atoi(optarg));
				break;

			case 'r':
				convolver_prefs.reuse_krn_fft = true;
				break;
#endif /* PROFIT_FFTW */

			case 'p':
				parse_profile(m, optarg);
				break;

#ifdef PROFIT_OPENCL
			case 'c':
				print_opencl_info();
				return 0;

			case 'C':
				use_opencl = true;
				tokenize(optarg, tokens, ",");
				if( tokens.size() != 3 ) {
					throw invalid_cmdline("-C argument must be of the form 'p,d,D' (e.g., -C 0,1,0)");
				}
				clplat_idx = (unsigned int)atoi(tokens[0].c_str());
				cldev_idx = (unsigned int)atoi(tokens[1].c_str());
				use_double = (bool)atoi(tokens[2].c_str());
				break;
#endif /* PROFIT_OPENCL */

#ifdef PROFIT_OPENMP
			case 'n':
				m.omp_threads = (unsigned int)atoi(optarg);
				convolver_prefs.omp_threads = m.omp_threads;
				break;
#endif

			case 'P':
				if( !stat(optarg, &stat_buf) ) {
					m.psf = read_image_from_fits_file(optarg, m.psf_width, m.psf_height, m.psf_scale_x, m.psf_scale_y);
				}
				else {
					m.psf = parse_psf(optarg, m.psf_width, m.psf_height, m.psf_scale_x, m.psf_scale_y);
				}
				convolver_prefs.krn_width = m.psf_width;
				convolver_prefs.krn_height = m.psf_height;
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

	m.width   = convolver_prefs.src_width = width;
	m.height  = convolver_prefs.src_height = height;
	m.scale_x = scale_x;
	m.scale_y = scale_y;
	m.magzero = magzero;

#ifdef PROFIT_OPENCL
	/* Get an OpenCL environment */
	if( use_opencl ) {
		auto start = system_clock::now();
		auto opencl_env = get_opencl_environment(clplat_idx, cldev_idx, use_double, show_stats);
		auto end = system_clock::now();
		m.opencl_env = opencl_env;
		convolver_prefs.opencl_env = opencl_env;
		auto opencl_duration = chrono::duration_cast<chrono::milliseconds>(end-start).count();
		cout << "OpenCL environment created in " << opencl_duration << " [ms]" << endl;
	}
#endif /* PROFIT_OPENCL */

	// Create the convolver
	auto start = system_clock::now();
	m.convolver = create_convolver(convolver_type, convolver_prefs);
	auto end = system_clock::now();

	auto duration = chrono::duration_cast<chrono::milliseconds>(end-start).count();
	cout << std::fixed << std::setprecision(3);
	cout << "Created convolver in " << duration << " [ms]" << endl;

	// Now run the model as many times as requested
	Image image = run(iterations, m);
	auto &imdata = image.getData();

	switch(output) {

		case none:
			break;

		case binary:
			fwrite(imdata.data(), sizeof(double), image.size(), stdout);
			break;

		case text:
			for(j=0; j!=image.getHeight(); j++) {
				for(i=0; i!=image.getWidth(); i++) {
					cout << imdata[j*image.getWidth() + i] << " ";
				}
				cout << endl;
			}
			break;

		case fits:
			if( to_fits(m, image, fits_output) ) {
				perror("Error while saving image to FITS file");
				return 1;
			}
			break;

		default:
			cerr << "Output not currently supported: " << output << endl;
	}

	if( show_stats ) {
		print_stats(m);
	}

	return 0;

}

int main(int argc, char *argv[]) {

#ifdef PROFIT_FFTW
	FFTPlan::initialize();
#endif

	int ret;
	try {
		ret = parse_and_run(argc, argv);
	}
	catch (invalid_cmdline &e) {
		cerr << "Error on command line: " << e.what() << endl;
		ret = 1;
	}
	catch (invalid_parameter &e) {
		cerr << "Error while calculating model: " << e.what() << endl;
		ret = 1;
	}
#ifdef PROFIT_OPENCL
	catch (opencl_error &e) {
		cerr << "Error in OpenCL operation: " << e.what() << endl;
		ret = 1;
	}
#endif /* PROFIT_OPENCL */
#ifdef PROFIT_FFTW
	catch (fft_error &e) {
		cerr << "Error in FFT operation: " << e.what() << endl;
		ret = 1;
	}
#endif /* PROFIT_FFT */

#ifdef PROFIT_FFTW
	FFTPlan::initialize();
#endif

	return ret;
}