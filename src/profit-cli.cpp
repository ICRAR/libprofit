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

// prevent min/max macros defined in windows.h to be defined in the first place
#ifdef _WIN32
# define NOMINMAX
# include <windows.h>
#endif // _WIN32

#include <getopt.h>

#include <algorithm>
#include <cstring>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <sstream>
#include <utility>

#include "profit/profit.h"
#include "profit/fits_utils.h"

namespace profit {

class invalid_cmdline : public std::exception {
public:
	explicit invalid_cmdline(const std::string& what) : m_what(what) {}
	invalid_cmdline(const invalid_cmdline &e) : m_what(e.m_what) {}
	~invalid_cmdline() throw() {}
	const char *what() const throw() { return m_what.c_str(); }

private:
	std::string m_what;
};

static
void parse_profile(std::ostream &os, Model &model, const std::string &description)
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
		try {
			p->parameter(parameter_spec);
		} catch (const unknown_parameter &e) {
			os << e.what();
		}
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
	unsigned int psf_width = stoui(*it++);
	unsigned int psf_height = stoui(*it++);
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

	auto stod = [] (const std::string &s) { return std::stod(s); };
	std::transform(values.begin(), values.end(), psf.begin(), stod);
	return psf;
}

static
void show_version(std::ostream &os) {
	using std::endl;
	os << "libprofit version " << version() << endl;
	os << "OpenCL support: ";
	if (has_opencl()) {
		os << "Yes (up to " << opencl_version_major() << "." << opencl_version_minor() << ")" << endl;
	}
	else {
		os << "No" << endl;
	}
	os << "OpenMP support: " << (has_openmp() ? "Yes" : "No") << endl;
	os << "FFTW support: ";
	if (has_fftw()) {
		os << "Yes ";
		if (has_fftw_with_openmp()) {
			os << "(with OpenMP)";
		}
		else {
			os << "(without OpenMP)";
		}
	}
	else {
		os << "No";
	}
	os << endl << "Extended CPU instruction sets supported:";
	bool sse2 = has_simd_instruction_set(simd_instruction_set::SSE2);
	bool avx = has_simd_instruction_set(simd_instruction_set::AVX);
	if (!sse2 && !avx) {
		os << " none";
	}
	if (sse2) {
		os << " SSE2";
	}
	if (avx) {
		os << " AVX";
	}
	os << endl;
}

static const char *help_msg = R"===(
%s: utility program to generate an image out of a model and a set of profiles

This program is licensed under the GPLv3 license.

Usage: %s [options] -p <spec> [-p <spec> ...]

Options:
  -t        Output image as text values on stdout
  -f <file> Output image as fits file
  -i <n>    Output performance information after evaluating the model n times
  -s        Show runtime stats
  -T <conv> Use this type of convolver (see below)
  -u        Return an un-cropped image from the convolver
  -C <p,d>  Use OpenCL with platform p, device d, and double support (0|1)
  -c        Display OpenCL information about devices and platforms
  -n <n>    Use n OpenMP threads to calculate profiles
  -e <n>    FFTW plans created with n effort (more takes longer)
  -I <n>    SIMD Instruction set to use with brute-force convolver.
            0=auto (default), 1=none, 2=sse2, 3=avx.
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

template <typename T>
static
void usage(std::basic_ostream<T> &os, char *prog_name) {
	char *buff = new char[std::strlen(help_msg) - 4 + std::strlen(prog_name) * 2 + 1];
	std::sprintf(buff, help_msg, prog_name, prog_name);
	os << buff;
	delete []buff;
}

static
void print_stats_line(std::ostream &os, const std::string &prefix, const std::string &stat_name, nsecs_t nsecs) {
	int nchars = int(prefix.size() + stat_name.size());
	int nspaces = std::max(0, 50 - nchars);
	std::string spaces(nspaces, ' ');
	os << prefix << stat_name << spaces << " : " << std::setw(10)
	   << std::setprecision(3) << std::fixed << double(nsecs) / 1e6 << " [ms]" << std::endl;
}

struct clver {
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
void print_opencl_info(std::ostream &out) {

	using std::endl;

	const auto info = get_opencl_info();

	if( info.size() > 0 ) {
		out << "OpenCL information" << endl;
		out << "==================" << endl << endl;
		for(auto platform_info: info) {
			auto plat_id = std::get<0>(platform_info);
			auto plat_info = std::get<1>(platform_info);
			out << "Platform [" << plat_id << "]" << endl;
			out << "  Name           : " << plat_info.name << endl;
			out << "  OpenCL version : " << clver{plat_info.supported_opencl_version} << endl;
			for(auto device_info: plat_info.dev_info) {
				out << "  Device [" << std::get<0>(device_info) << "]" << endl;
				out << "    Name           : " << std::get<1>(device_info).name << endl;
				out << "    OpenCL version : " << clver {std::get<1>(device_info).cl_version} << endl;
				out << "    Double         : " << (std::get<1>(device_info).double_support ? "Supported" : "Not supported") << endl;
			}
			out << endl;
		}
	}
	else {
		out << "No OpenCL installation found" << endl;
	}
}

static
void print_cl_command_times(std::ostream &os, const std::string &prefix, const OpenCL_command_times &t, const std::string &action)
{
	print_stats_line(os, prefix, action + " submission", t.submit);
	print_stats_line(os, prefix, action + " waiting", t.wait);
	print_stats_line(os, prefix, action + " execution", t.exec);
}

static
void print_cl_stats(std::ostream &os, const std::string &prefix0, bool opencl_120, const OpenCL_times &stats) {

	auto prefix1 = prefix0 + "  ";

	std::ostringstream cl_ops_os;
	cl_ops_os << "OpenCL operations (" << stats.nwork_items << " work items)";
	print_stats_line(os, prefix0, cl_ops_os.str(), stats.total);
	print_stats_line(os, prefix1, "Kernel preparation", stats.kernel_prep);
	if( opencl_120 ) {
		print_cl_command_times(os, prefix1, stats.filling_times, "Fill");
	}
	print_cl_command_times(os, prefix1, stats.writing_times, "Write");
	print_cl_command_times(os, prefix1, stats.kernel_times, "Kernel");
	print_cl_command_times(os, prefix1, stats.reading_times, "Read");
}

static
void print_stats(std::ostream &os, const Model &m) {

#ifdef PROFIT_DEBUG
	for(const auto &profile_integrations: m.get_profile_integrations()) {
		int total = 0;
		if( std::get<1>(profile_integrations).size() > 0 ) {
			os << "Integrations per recursion level for profile " << std::get<0>(profile_integrations) << std::endl;
			for(const auto level_integrations: std::get<1>(profile_integrations)) {
				auto integrations = std::get<1>(level_integrations);
				total += integrations;
				os << " Level " << std::get<0>(level_integrations) << ": " << integrations << " integrations" << std::endl;
			}
			os << " Total: " << total << " integrations" << std::endl;
		}
		else {
			os << "Profile " << std::get<0>(profile_integrations) << " didn't run into any recursion" << std::endl;
		}
	}
#endif /* PROFIT_DEBUG */

	os << std::endl;
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

		os << "Stats for profile " << profile_name << std::endl;

		auto prefix1 = "  ";
		RadialProfileStats *rprofile_stats = dynamic_cast<RadialProfileStats *>(profile_stats);
		auto opencl_env = m.get_opencl_env();
		if( rprofile_stats && opencl_env ) {
			bool opencl_120 = opencl_env->get_version() >= 120;
			print_cl_stats(os, prefix0, opencl_120, rprofile_stats->cl_times);
			print_stats_line(os, prefix0, "Pre-loop", rprofile_stats->subsampling.pre_subsampling);
			print_stats_line(os, prefix0, "Subsampling loop", rprofile_stats->subsampling.total);
			print_stats_line(os, prefix1, "New subsamples calculation", rprofile_stats->subsampling.new_subsampling);
			print_stats_line(os, prefix1, "Initial transform", rprofile_stats->subsampling.inital_transform);
			print_cl_stats(os, prefix1, opencl_120, rprofile_stats->subsampling.cl_times);
			print_stats_line(os, prefix1, "Final transform", rprofile_stats->subsampling.final_transform);
			print_stats_line(os, prefix0, "Final image", rprofile_stats->final_image);
		}

		print_stats_line(os, prefix0, "Total", profile_stats->total);
	}
}

static
Image run(std::ostream &os, unsigned int iterations, Model &m, Point &offset) {

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
	os << std::fixed << std::setprecision(3);
	os << "Ran " << iterations << " iterations in ";
	os << std::setprecision(3) << std::fixed << dur_secs << " [s] ";
	os << "(" << std::setprecision(3) << std::fixed << dur_per_iter << " [ms] per iteration)";
	os << std::endl;

	return result;
}

typedef enum _output_type {
	none = 0,
	text,
	fits,
} output_t;

static
int parse_and_run(int argc, char *argv[], std::ostream &cout, std::ostream &cerr) {

	namespace chrono = std::chrono;
	using chrono::system_clock;

	int opt;
	unsigned int width = 100;
	unsigned int height = 100;
	unsigned int iterations = 1;
	double scale_x = 1;
	double scale_y = 1;
	unsigned int i;
	unsigned int j;
	std::string fits_output;
	output_t output = none;
	Model m;
	Image psf;
	unsigned int finesampling = 1;
	std::string convolver_type = "brute";
	ConvolverCreationPreferences convolver_prefs;
	bool show_stats = false;

	bool use_opencl = false;
	bool use_double = false;
	unsigned int clplat_idx = 0;
	unsigned int cldev_idx = 0;
	std::vector<std::string> tokens;

	const char *options = "h?VsRP:p:w:H:x:y:X:Y:m:tf:i:T:uS:C:ce:rn:FI:";

	while( (opt = getopt(argc, argv, options)) != -1 ) {
		switch(opt) {

			case 'h':
			case '?':
				usage(cout, argv[0]);
				return 0;

			case 'V':
				show_version(cout);
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
				convolver_prefs.effort = effort_t(std::stoul(optarg));
				break;

			case 'r':
				convolver_prefs.reuse_krn_fft = true;
				break;

			case 'I':
				convolver_prefs.instruction_set = simd_instruction_set(std::stoul(optarg));
				break;

			case 'p':
				parse_profile(cerr, m, optarg);
				break;

			case 'c':
				print_opencl_info(cout);
				return 0;

			case 'C':
				if (!has_opencl()) {
					throw invalid_cmdline("libprofit was compiled without OpenCL support, but support was requested. See -V for details");
				}
				use_opencl = true;
				tokens = split(optarg, ",");
				if( tokens.size() != 3 ) {
					throw invalid_cmdline("-C argument must be of the form 'p,d,D' (e.g., -C 0,1,0)");
				}
				clplat_idx = stoui(tokens[0].c_str());
				cldev_idx = stoui(tokens[1].c_str());
				use_double = bool(stoui(tokens[2].c_str()));
				break;

			case 'n':
				m.set_omp_threads(stoui(optarg));
				convolver_prefs.omp_threads = m.get_omp_threads();
				break;

			case 'P':
				if( file_exists(optarg) ) {
					try {
						PixelScale psf_pixel_scale;
						psf = from_fits(optarg, psf_pixel_scale);
						m.set_psf_pixel_scale(psf_pixel_scale);
					} catch (const invalid_file &e) {
						std::ostringstream os;
						os << "Error loading PSF from " << optarg << ": " << e.what();
						throw invalid_cmdline(os.str());
					}
				}
				else {
					psf = parse_psf(optarg, m);
				}
				convolver_prefs.krn_dims = psf.getDimensions();
				m.set_psf(std::move(psf));
				break;

			case 'w':
				width = stoui(optarg);
				break;

			case 'H':
				height = stoui(optarg);
				break;

			case 'S':
				finesampling = stoui(optarg);
				break;

			case 'F':
				m.set_return_finesampled(false);
				break;

			case 'x':
				scale_x = std::stod(optarg);
				break;

			case 'y':
				scale_y = std::stod(optarg);
				break;

			case 'm':
				m.set_magzero(std::stod(optarg));
				break;

			case 't':
				output = text;
				break;

			case 'f':
				fits_output = optarg;
				output = fits;
				break;

			case 'i':
				iterations = stoui(optarg);
				break;

			default:
				usage(cerr, argv[0]);
				return 1;

		}
	}

	/* No profiles given */
	if( !m.has_profiles() ) {
		usage(cerr, argv[0]);
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
		cout << "OpenCL environment (platform=" <<
		        opencl_env->get_platform_name() << ", device=" <<
		        opencl_env->get_device_name() << ", version=" <<
		        clver{opencl_env->get_version()} <<
		        ") created in " << opencl_duration << " [ms]" << std::endl;
	}

	// Create the convolver
	auto start = system_clock::now();
	m.set_convolver(create_convolver(convolver_type, convolver_prefs));
	auto end = system_clock::now();

	auto duration = chrono::duration_cast<chrono::milliseconds>(end-start).count();
	cout << std::fixed << std::setprecision(3);
	cout << "Created convolver in " << duration << " [ms]" << std::endl;

	// Now run the model as many times as requested
	Point offset;
	auto image = run(cout, iterations, m, offset);

	switch(output) {

		case none:
			break;

		case text:
			for(j=0; j!=image.getHeight(); j++) {
				for(i=0; i!=image.getWidth(); i++) {
					cout << image[j*image.getWidth() + i] << " ";
				}
				cout << std::endl;
			}
			break;

		case fits:
			to_fits(image, offset, m.get_image_pixel_scale(), fits_output);
			break;

		default:
			cerr << "Output not currently supported: " << output << std::endl;
	}

	if( show_stats ) {
		print_stats(cout, m);
	}

	return 0;
}

} // namespace shark

extern "C" {

int main(int argc, char *argv[]) {

	std::ostream &cout = std::cout;
	std::ostream &cerr = std::cerr;

	bool success = profit::init();
	auto init_diagnose = profit::init_diagnose();
	if (!success) {
		cerr << "Error initializing libprofit: " << init_diagnose << std::endl;
		return 1;
	}
	else if (!init_diagnose.empty()){
		cerr << "Warning while initializing libprofit: " << init_diagnose << std::endl;
	}

	int ret;
	try {
		ret = profit::parse_and_run(argc, argv, cout, cerr);
	}
	catch (const profit::invalid_cmdline &e) {
		cerr << "Error on command line: " << e.what() << std::endl;
		ret = 1;
	}
	catch (const profit::invalid_parameter &e) {
		cerr << "Error while calculating model: " << e.what() << std::endl;
		ret = 1;
	}
	catch (const profit::opencl_error &e) {
		cerr << "Error in OpenCL operation: " << e.what() << std::endl;
		ret = 1;
	}
	catch (const profit::fft_error &e) {
		cerr << "Error in FFT operation: " << e.what() << std::endl;
		ret = 1;
	}
	catch (const std::exception &e) {
		cerr << "Unexpected error: " << e.what() << std::endl;
		ret = 1;
	}
	profit::finish();
	auto finish_diagnose = profit::finish_diagnose();
	if (!finish_diagnose.empty()) {
		cerr << "Warning while finishing libprofit: " << finish_diagnose << std::endl;
	}

	return ret;
}

} // extern "C"