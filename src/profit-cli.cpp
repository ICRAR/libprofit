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
#include <iterator>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <sstream>

#include "profit/profit.h"


namespace profit {

class invalid_cmdline : std::exception {
public:
	invalid_cmdline(const std::string& what) : m_what(what) {}
	invalid_cmdline(const invalid_cmdline &e) : m_what(e.m_what) {}
	~invalid_cmdline() throw() {}
	const char *what() const throw() { return m_what.c_str(); }

private:
	std::string m_what;
};

static
void parse_profile(Model &model, const std::string &description)
{
	auto desc = trim(description);
	if (desc.empty()) {
		throw invalid_cmdline("Missing parameter name after -p");
	}

	/* The description might be only a name */
	auto parts = split(desc, ":");
	auto p = model.add_profile(parts[0]);
	std::vector<std::string> parameter_specs(std::make_move_iterator(parts.begin() + 1),
	                                         std::make_move_iterator(parts.end()));
	for(auto &parameter_spec: parameter_specs) {
		p->parameter(parameter_spec);
	}
}

static
Image parse_psf(std::string optarg, Model &m)
{
	bool read_scales = false;

	/* format is w:h:[optional scale_x:scale_y:]:val1,val2... */
	auto tokens = split(optarg, ":");
	auto ntokens = tokens.size();
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

	auto it = tokens.begin();
	unsigned int psf_width = std::stoul(*it++);
	unsigned int psf_height = std::stoul(*it++);
	if( read_scales ) {
		double psf_scale_x = std::stod(*it++);
		double psf_scale_y = std::stod(*it++);
		m.set_psf_pixel_scale({psf_scale_x, psf_scale_y});
	}

	Image psf(psf_width, psf_height);
	auto values = split(*it, ",");
	if (values.size() != psf.size()) {
		std::ostringstream os;
		os << "Not enough values provided for PSF. Provided: " << values.size() << ", expected: " << psf.size();
		throw invalid_cmdline(os.str());
	}

	auto stod = [] (const std::string &s) -> double { return std::stod(s); };
	std::transform(values.begin(), values.end(), psf.begin(), stod);
	return psf;
}

static
void show_version() {
	using std::cout;
	using std::endl;
	cout << "libprofit version " << version() << endl;
	cout << "OpenCL support: ";
	if (has_opencl()) {
		cout << "Yes (up to " << opencl_version_major() << "." << opencl_version_minor() << ")" << endl;
	}
	else {
		cout << "No" << endl;
	}
	cout << "OpenMP support: " << (has_openmp() ? "Yes" : "No") << endl;
	cout << "FFTW support: ";
	if (has_fftw()) {
		cout << "Yes ";
		if (has_fftw_with_openmp()) {
			cout << "(with OpenMP)";
		}
		else {
			cout << "(without OpenMP)";
		}
	}
	else {
		cout << "No";
	}
	cout << endl;
}

static const char *help_msg = R"===(
%s: utility program to generate an image out of a model and a set of profiles

This program is licensed under the GPLv3 license.

Usage: %s [options] -p <spec> [-p <spec> ...]

Options:
  -t        Output image as text values on stdout
  -b        Output image as binary content on stdout
  -f <file> Output image as fits file
  -i <n>    Output performance information after evaluating the model n times
  -s        Show runtime stats
  -T <conv> Use this type of convolver (see below)
  -u        Return an un-cropped image from the convolver
  -C <p,d>  Use OpenCL with platform p, device d, and double support (0|1)
  -c        Display OpenCL information about devices and platforms
  -n <n>    Use n OpenMP threads to calculate profiles
  -e <n>    FFTW plans created with n effort (more takes longer)
  -r        Reuse FFT-transformed PSF across evaluations (if -T fft)
  -x        Image width. Defaults to 100
  -y        Image height. Defaults to 100
  -S <n>    Finesampling factor. Defaults to 1
  -F        Do *not* return finesampled image (if -S <n>)
  -w        Width in pixels. Defaults to 100
  -H        Height in pixels. Defaults to 100
  -m        Zero magnitude. Defaults to 0
  -P        PSF function (specified as w:h:val1,val2..., or as a FITS filename)
  -R        Clear libprofit's cache and exit
  -h,-?     Show this help and exit
  -V        Show the program version and exit

The following convolver types are supported:

 * brute: A brute-force convolver
 * brute-old: An older, slower brute-force convolver (used only for comparisons)
 * opencl: An OpenCL-based brute-force convolver
 * fft: An FFT-based convolver

Profiles should be specified as follows:

-p name:param1=val1:param2=val2:...

The following profiles (and parameters) are currently accepted:

 * psf: xcen, ycen, mag
 * sky: bg
 * sersic: re, nser, rescale_flux
 * moffat: fwhm, con
 * ferrer: a, b, rout
 * coresersic: re, nser, rb, a, b
 * brokenexp: h1, h2, rb, a
 * king: rc, rt, a
 * sersic, moffat, ferrer, coresersic, king: xcen, ycen, mag, box, ang, axrat,
                           rough, rscale_switch, max_recursions,
                           resolution, acc, rscale_max, adjust

For more information visit https://libprofit.readthedocs.io.

)===";

static
void usage(FILE *file, char *argv[]) {
	std::fprintf(file, help_msg, argv[0], argv[0]);
}

static
void print_stats_line(const std::string &prefix, const std::string &stat_name, double val) {
	const auto static name_width = 50u;
	int nchars = prefix.size() + stat_name.size();
	int nspaces = name_width - nchars;
	std::string spaces(std::max(0, nspaces), ' ');
	std::cout << prefix << stat_name << spaces << " : " << std::setw(10)
	          << std::setprecision(3) << std::fixed << val << " [ms]" << std::endl;
}

struct clver {
	clver(unsigned int ver) : ver(ver) {}
	unsigned int ver;
};

template <typename T>
static
std::basic_ostream<T> &operator<<(std::basic_ostream<T> &os, const clver &ver)
{
	auto major = ver.ver / 100;
	auto minor = (ver.ver - major * 100) / 10;
	os << major << "." << minor;
	return os;
}

static
void print_opencl_info() {

	using std::cout;
	using std::endl;

	const auto info = get_opencl_info();

	if( info.size() > 0 ) {
		cout << "OpenCL information" << endl;
		cout << "==================" << endl << endl;
		for(auto platform_info: info) {
			auto plat_id = std::get<0>(platform_info);
			auto plat_info = std::get<1>(platform_info);
			cout << "Platform [" << plat_id << "]" << endl;
			cout << "  Name           : " << plat_info.name << endl;
			cout << "  OpenCL version : " << clver(plat_info.supported_opencl_version) << endl;
			for(auto device_info: plat_info.dev_info) {
				cout << "  Device [" << std::get<0>(device_info) << "]" << endl;
				cout << "    Name           : " << std::get<1>(device_info).name << endl;
				cout << "    OpenCL version : " << clver(std::get<1>(device_info).cl_version) << endl;
				cout << "    Double         : " << (std::get<1>(device_info).double_support ? "Supported" : "Not supported") << endl;
			}
			cout << endl;
		}
	}
	else {
		cout << "No OpenCL installation found" << endl;
	}
}

static
void print_cl_stats(const std::string &prefix0, bool opencl_120, const OpenCL_times &stats) {

	auto prefix1 = prefix0 + "  ";

	std::ostringstream os;
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

static
void print_stats(const Model &m) {

#ifdef PROFIT_DEBUG
	for(const auto &profile_integrations: m.get_profile_integrations()) {
		int total = 0;
		if( std::get<1>(profile_integrations).size() > 0 ) {
			std::cout << "Integrations per recursion level for profile " << std::get<0>(profile_integrations) << std::endl;
			for(const auto level_integrations: std::get<1>(profile_integrations)) {
				auto integrations = std::get<1>(level_integrations);
				total += integrations;
				std::cout << " Level " << std::get<0>(level_integrations) << ": " << integrations << " integrations" << std::endl;
			}
			std::cout << " Total: " << total << " integrations" << std::endl;
		}
		else {
			std::cout << "Profile " << std::get<0>(profile_integrations) << " didn't run into any recursion" << std::endl;
		}
	}
#endif /* PROFIT_DEBUG */

	std::cout << std::endl;
	auto const &stats = m.get_stats();

	auto prefix0 = "";
	for(auto const &stat_pair: stats) {

		// Some profile might not have gathered stats
		auto profile_name = std::get<0>(stat_pair);
		auto stat = std::get<1>(stat_pair);
		ProfileStats *profile_stats = stat.get();
		if( !profile_stats ) {
			continue;
		}

		std::cout << "Stats for profile " << profile_name << std::endl;

		auto prefix1 = "  ";
		RadialProfileStats *rprofile_stats = dynamic_cast<RadialProfileStats *>(profile_stats);
		auto opencl_env = m.get_opencl_env();
		if( rprofile_stats && opencl_env ) {
			bool opencl_120 = opencl_env->get_version() >= 120;
			print_cl_stats(prefix0, opencl_120, rprofile_stats->cl_times);
			print_stats_line(prefix0, "Pre-loop", rprofile_stats->subsampling.pre_subsampling / 1e6 );
			print_stats_line(prefix0, "Subsampling loop", rprofile_stats->subsampling.total / 1e6 );
			print_stats_line(prefix1, "New subsamples calculation", rprofile_stats->subsampling.new_subsampling / 1e6 );
			print_stats_line(prefix1, "Initial transform", rprofile_stats->subsampling.inital_transform / 1e6 );
			print_cl_stats(prefix1, opencl_120, rprofile_stats->subsampling.cl_times);
			print_stats_line(prefix1, "Final transform", rprofile_stats->subsampling.final_transform / 1e6 );
			print_stats_line(prefix0, "Final image", rprofile_stats->final_image / 1e6 );
		}

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

static
Image read_image_from_fits_file(const std::string &filename, Model &m) {

	FILE *f;
	unsigned int pos, padding;
	char hdr[80];

	unsigned int width = 0;
	unsigned int height = 0;
	double scale_x = 1;
	double scale_y = 1;

	f = fopen(filename.c_str(), "rb");
	if( !f ) {
		std::ostringstream os;
		os << "Couldn't open '" << filename << "' for reading: " << strerror(errno);
		throw invalid_cmdline(os.str());
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

	m.set_psf_pixel_scale({scale_x, scale_y});

	pos = (unsigned int)ftell(f);
	padding = FITS_BLOCK_SIZE - (pos % FITS_BLOCK_SIZE);
	fseek(f, padding, SEEK_CUR);

	Image psf(width, height);
	size_t nitems = fread(psf.data(), sizeof(double), psf.size(), f);
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

int to_fits(Model &m, const Image &image, const Point &offset, std::string fname) {

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
	auto scale_x = m.get_image_pixel_scale().first;
	auto scale_y = m.get_image_pixel_scale().second;
	fprintf(f, "%-80s", "SIMPLE  =                    T / File conforms to FITS standard");
	fprintf(f, "%-80s", "BITPIX  =                  -64 / Bits per pixel");
	fprintf(f, "%-80s", "NAXIS   =                    2 / Number of axes");
	sprintf(hdr, "NAXIS1  =           %10.0u / Width", image.getWidth());
	fprintf(f, "%-80s", hdr);
	sprintf(hdr, "NAXIS2  =           %10.0u / Height", image.getHeight());
	fprintf(f, "%-80s", hdr);
	fprintf(f, "%-80s", "CRPIX1  = 1");
	sprintf(hdr, "CRVAL1  = %f", (0.5 - offset.x) * scale_x);
	fprintf(f, "%-80s", hdr);
	sprintf(hdr, "CDELT1  = %f", scale_x);
	fprintf(f, "%-80s", hdr);
	fprintf(f, "%-80s", "CTYPE1  = ' '");
	fprintf(f, "%-80s", "CUNIT1  = ' '");
	fprintf(f, "%-80s", "CRPIX2  = 1");
	sprintf(hdr, "CRVAL2  = %f", (0.5 - offset.y) * scale_y);
	fprintf(f, "%-80s", hdr);
	sprintf(hdr, "CDELT2  = %f", scale_y);
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
	if( is_little_endian() ) {
		std::vector<double> big_endian_image(image_size);
		std::transform(image.begin(), image.end(), big_endian_image.begin(), swap_bytes);
		fwrite(big_endian_image.data(), sizeof(double), image_size, f);
	}
	else {
		fwrite(image.data(), sizeof(double), image_size, f);
	}

	/* Pad with zeroes until we complete the current 36*80 block */
	padding = FITS_BLOCK_SIZE - (((unsigned int)sizeof(double) * image_size) % FITS_BLOCK_SIZE);
	std::string zeros(padding, 0);
	fwrite(zeros.c_str(), 1, padding, f);
	fclose(f);

	return 0;
}

static
Image run(unsigned int iterations, Model &m, Point &offset) {

	using std::chrono::system_clock;

	/* This means that we evaluated the model once, but who cares */
	Image result;
	auto start = system_clock::now();
	for(unsigned i=0; i!=iterations; i++) {
		result = m.evaluate(offset);
	}
	auto end = system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

	double dur_secs = (double)duration/1000;
	double dur_per_iter = (double)duration/iterations;
	std::cout << std::fixed << std::setprecision(3);
	std::cout << "Ran " << iterations << " iterations in ";
	std::cout << std::setprecision(3) << std::fixed << dur_secs << " [s] ";
	std::cout << "(" << std::setprecision(3) << std::fixed << dur_per_iter << " [ms] per iteration)";
	std::cout << std::endl;

	return result;
}

typedef enum _output_type {
	none = 0,
	binary = 1,
	text = 2,
	fits = 3,
} output_t;

static
int parse_and_run(int argc, char *argv[]) {

	namespace chrono = std::chrono;
	using chrono::system_clock;

	int opt;
	unsigned int width = 100, height = 100, iterations = 1;
	double scale_x = 1, scale_y = 1;
	unsigned int i, j;
	char *endptr = NULL;
	std::string fits_output;
	output_t output = none;
	Model m;
	Image psf;
	unsigned int finesampling = 1;
	std::string convolver_type = "brute";
	ConvolverCreationPreferences convolver_prefs;
	struct stat stat_buf;
	bool show_stats = false;

	bool use_opencl = false, use_double = false;
	unsigned int clplat_idx = 0, cldev_idx = 0;
	std::vector<std::string> tokens;

	const char *options = "h?VsRP:p:w:H:x:y:X:Y:m:tbf:i:T:uS:C:ce:rn:F";

	while( (opt = getopt(argc, argv, options)) != -1 ) {
		switch(opt) {

			case 'h':
			case '?':
				usage(stdout, argv);
				return 0;

			case 'V':
				show_version();
				return 0;

			case 'R':
				clear_cache();
				return 0;

			case 's':
				show_stats = true;
				break;

			case 'T':
				convolver_type = optarg;
				break;

			case 'u':
				m.set_crop(false);
				break;

			case 'e':
				convolver_prefs.effort = effort_t(std::atoi(optarg));
				break;

			case 'r':
				convolver_prefs.reuse_krn_fft = true;
				break;

			case 'p':
				parse_profile(m, optarg);
				break;

			case 'c':
				print_opencl_info();
				return 0;

			case 'C':
				if (not has_opencl()) {
					throw invalid_cmdline("libprofit was compiled without OpenCL support, but support was requested. See -V for details");
				}
				use_opencl = true;
				tokens = split(optarg, ",");
				if( tokens.size() != 3 ) {
					throw invalid_cmdline("-C argument must be of the form 'p,d,D' (e.g., -C 0,1,0)");
				}
				clplat_idx = (unsigned int)atoi(tokens[0].c_str());
				cldev_idx = (unsigned int)atoi(tokens[1].c_str());
				use_double = (bool)atoi(tokens[2].c_str());
				break;

			case 'n':
				m.set_omp_threads((unsigned int)atoi(optarg));
				convolver_prefs.omp_threads = m.get_omp_threads();
				break;

			case 'P':
				if( !stat(optarg, &stat_buf) ) {
					psf = read_image_from_fits_file(optarg, m);
				}
				else {
					psf = parse_psf(optarg, m);
				}
				convolver_prefs.krn_dims = psf.getDimensions();
				m.set_psf(std::move(psf));
				break;

			case 'w':
				width = (unsigned int)atoi(optarg);
				break;

			case 'H':
				height = (unsigned int)atoi(optarg);
				break;

			case 'S':
				finesampling = (unsigned int)atoi(optarg);
				break;

			case 'F':
				m.set_return_finesampled(false);
				break;

			case 'x':
				scale_x = atof(optarg);
				break;

			case 'y':
				scale_y = atof(optarg);
				break;

			case 'm':
				m.set_magzero(strtod(optarg, &endptr));
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

	Dimensions dims {width, height};
	m.set_dimensions(dims);
	m.set_image_pixel_scale({scale_x, scale_y});
	m.set_finesampling(finesampling);
	convolver_prefs.src_dims = dims * finesampling;

	/* Get an OpenCL environment */
	if( use_opencl ) {
		auto start = system_clock::now();
		auto opencl_env = get_opencl_environment(clplat_idx, cldev_idx, use_double, show_stats);
		auto end = system_clock::now();
		m.set_opencl_env(opencl_env);
		convolver_prefs.opencl_env = opencl_env;
		auto opencl_duration = chrono::duration_cast<chrono::milliseconds>(end-start).count();
		std::cout << "OpenCL environment (platform=" <<
		            opencl_env->get_platform_name() << ", device=" <<
		            opencl_env->get_device_name() << ", version=" <<
		            clver(opencl_env->get_version()) <<
		            ") created in " << opencl_duration << " [ms]" << std::endl;
	}

	// Create the convolver
	auto start = system_clock::now();
	m.set_convolver(create_convolver(convolver_type, convolver_prefs));
	auto end = system_clock::now();

	auto duration = chrono::duration_cast<chrono::milliseconds>(end-start).count();
	std::cout << std::fixed << std::setprecision(3);
	std::cout << "Created convolver in " << duration << " [ms]" << std::endl;

	// Now run the model as many times as requested
	Point offset;
	auto image = run(iterations, m, offset);

	switch(output) {

		case none:
			break;

		case binary:
			fwrite(image.data(), sizeof(double), image.size(), stdout);
			break;

		case text:
			for(j=0; j!=image.getHeight(); j++) {
				for(i=0; i!=image.getWidth(); i++) {
					std::cout << image[j*image.getWidth() + i] << " ";
				}
				std::cout << std::endl;
			}
			break;

		case fits:
			if( to_fits(m, image, offset, fits_output) ) {
				perror("Error while saving image to FITS file");
				return 1;
			}
			break;

		default:
			std::cerr << "Output not currently supported: " << output << std::endl;
	}

	if( show_stats ) {
		print_stats(m);
	}

	return 0;
}

} // namespace shark

extern "C" {

int main(int argc, char *argv[]) {

	bool success = profit::init();
	auto init_diagnose = profit::init_diagnose();
	if (!success) {
		std::cerr << "Error initializing libprofit: " << init_diagnose << std::endl;
		return 1;
	}
	else if (!init_diagnose.empty()){
		std::cerr << "Warning while initializing libprofit: " << init_diagnose << std::endl;
	}

	int ret;
	try {
		ret = profit::parse_and_run(argc, argv);
	}
	catch (const profit::invalid_cmdline &e) {
		std::cerr << "Error on command line: " << e.what() << std::endl;
		ret = 1;
	}
	catch (const profit::invalid_parameter &e) {
		std::cerr << "Error while calculating model: " << e.what() << std::endl;
		ret = 1;
	}
	catch (const profit::opencl_error &e) {
		std::cerr << "Error in OpenCL operation: " << e.what() << std::endl;
		ret = 1;
	}
	catch (const profit::fft_error &e) {
		std::cerr << "Error in FFT operation: " << e.what() << std::endl;
		ret = 1;
	}
	catch (const std::exception &e) {
		std::cerr << "Unexpected error: " << e.what() << std::endl;
		ret = 1;
	}
	profit::finish();
	auto finish_diagnose = profit::finish_diagnose();
	if (!finish_diagnose.empty()) {
		std::cerr << "Warning while finishing libprofit: " << finish_diagnose << std::endl;
	}

	return ret;
}

} // extern "C"