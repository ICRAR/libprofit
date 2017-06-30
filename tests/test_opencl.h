/**
 * OpenCL-related tests
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2016
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Rodrigo Tobar
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

#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include <cxxtest/GlobalFixture.h>
#include <cxxtest/TestSuite.h>

#include "profit/profit.h"

#ifdef PROFIT_OPENCL

using namespace profit;

class OpenCLFixtures : CxxTest::GlobalFixture {

public:
	OpenCLFixtures() {

		int plat_idx = -1, dev_idx = -1;
		bool use_double;

		// Look preferably for an OpenCL device that has double support
		try {
			auto platforms = get_opencl_info();
			for(auto &platform: platforms) {
				plat_idx = std::get<0>(platform);
				auto platform_info = std::get<1>(platform);
				for(auto dev: platform_info.dev_info) {
					dev_idx = std::get<0>(dev);
					auto dev_info = std::get<1>(dev);
					use_double = dev_info.double_support;
					if ( use_double ) {
						goto chosen;
					}
				}
			}
		} catch (const opencl_error &e) {
			// no-op
		}

	chosen:
		// Actually nothing found, or there was an error while trying
		if ( plat_idx == -1 ) {
			return;
		}

		try {
			opencl_env = get_opencl_environment(plat_idx, dev_idx, use_double, false);
		} catch (const opencl_error &e) {
			// we simply loose hope
			return;
		}

		// 5% tolerance by default, can be overridden via environment variable
		tolerance = 0.05;
		if( const char *tol = std::getenv("LIBPROFIT_OPENCL_TOLERANCE") ) {
			tolerance = std::stod(tol);
		}

	}

	double tolerance;
	std::shared_ptr<OpenCL_env> opencl_env;
};
static OpenCLFixtures openCLFixtures;

class TestOpenCL : public CxxTest::TestSuite {

private:

	void _pixels_within_tolerance(std::vector<double> original_im, std::vector<double> opencl_im,
	                              unsigned int i, unsigned int width,
	                              double tolerance) {

		double original = original_im[i];
		double opencl = opencl_im[i];
		auto diff = std::abs(original - opencl);
		if ( !diff ) {
			// all good
			return;
		}

		// avoid NaNs due to divide-by-zero
		auto denomin = original;
		if ( !denomin ) {
			denomin = opencl;
		}
		auto relative_diff = diff / denomin;

		std::ostringstream msg;
		msg << "Pixel [" << i % width << "," << i / width << "] has values that are too different: ";
		msg << original << " v/s " << opencl;
		TSM_ASSERT_LESS_THAN_EQUALS(msg.str(), relative_diff, tolerance);
	}

	void _check_images_within_tolerance(Model &m) {

		if( !openCLFixtures.opencl_env ) {
			TS_SKIP("No OpenCL environment found to run OpenCL tests");
		}

		// evaluate normally first, and then using the OpenCL environment,
		// which they all support
		std::vector<double> original = m.evaluate();
		m.opencl_env = openCLFixtures.opencl_env;
		std::vector<double> opencl_produced = m.evaluate();

		// Pixel by pixel the images should be fairly similar
		for(unsigned int i=0; i!=original.size(); i++) {
			_pixels_within_tolerance(original, opencl_produced, i, m.width, openCLFixtures.tolerance);
		}
	}

	void _check_convolver(Convolver &&clConvolver) {
		for(auto im_dim: {100, 101}) {
			for(auto krn_dim: {24, 25}) {
				Mask mask;
				Image src(im_dim, im_dim);
				Image krn(krn_dim, krn_dim);
				for(auto &d: src.getData()) {
					d = (rand() % 10000) / 10000.0;
				}
				for(auto &d: krn.getData()) {
					d = (rand() % 10000) / 10000.0;
				}

				BruteForceConvolver bConvolver;
				Image result1 = bConvolver.convolve(src, krn, mask);
				Image result2 = clConvolver.convolve(src, krn, mask);
				for(unsigned int i = 0; i < src.getSize(); i++) {
					// Hopefully within 0.1% of error?
					_pixels_within_tolerance(result1.getData(), result2.getData(), i, src.getWidth(), 1e-3);
				}
			}
		}
	}

public:

	void test_opencldiff_brokenexp() {
		Model m;
		m.width = m.height = 100;
		auto brokenexp = m.add_profile("brokenexp");
		brokenexp->parameter("xcen", 40.);
		brokenexp->parameter("ycen", 60.);
		brokenexp->parameter("h1", 10.);
		brokenexp->parameter("h2", 5.0);
		brokenexp->parameter("rb", 8.0);
		brokenexp->parameter("ang", 30.);
		brokenexp->parameter("a", 0.2);
		brokenexp->parameter("axrat", 0.6);
		brokenexp->parameter("box", -0.1);
		_check_images_within_tolerance(m);
	}

	void test_opencldiff_coresersic() {
		Model m;
		m.width = m.height = 100;
		auto coresersic = m.add_profile("coresersic");
		coresersic->parameter("xcen", 50.);
		coresersic->parameter("ycen", 50.);
		coresersic->parameter("mag", 15.);
		coresersic->parameter("rb", 5.0);
		coresersic->parameter("re", 10.0);
		coresersic->parameter("nser", 4.0);
		coresersic->parameter("a", 2.);
		coresersic->parameter("b", 1.3);
		coresersic->parameter("ang", 30.);
		coresersic->parameter("axrat", 0.4);
		coresersic->parameter("box", 0.);
		_check_images_within_tolerance(m);
	}

	void test_opencldiff_ferrer() {
		Model m;
		m.width = m.height = 100;
		auto ferrer = m.add_profile("ferrer");
		ferrer->parameter("xcen", 50.);
		ferrer->parameter("ycen", 50.);
		ferrer->parameter("mag", 15.);
		ferrer->parameter("rout", 20.0);
		ferrer->parameter("a", 0.5);
		ferrer->parameter("b", 0.7);
		ferrer->parameter("ang", 30.);
		ferrer->parameter("axrat", 0.3);
		ferrer->parameter("box", 0.6);
		_check_images_within_tolerance(m);
	}

	void test_opencldiff_king() {
		Model m;
		m.width = m.height = 100;
		auto king = m.add_profile("king");
		king->parameter("xcen", 50.);
		king->parameter("ycen", 50.);
		king->parameter("mag", 15.);
		king->parameter("rc", 5.0);
		king->parameter("rt", 30.);
		king->parameter("a", 2.0);
		king->parameter("ang", 30.);
		king->parameter("axrat", 0.8);
		king->parameter("box", 0.0);
		_check_images_within_tolerance(m);
	}

	void test_opencldiff_moffat() {
		Model m;
		m.width = m.height = 100;
		auto moffat = m.add_profile("moffat");
		moffat->parameter("xcen", 34.);
		moffat->parameter("ycen", 74.);
		moffat->parameter("mag", 1.);
		moffat->parameter("fwhm", 3.0);
		moffat->parameter("con", 5.);
		moffat->parameter("ang", 45.);
		moffat->parameter("axrat", 0.95);
		moffat->parameter("box", 0.1);
		_check_images_within_tolerance(m);
	}

	void test_opencldiff_sersic() {
		Model m;
		m.width = m.height = 100;
		auto sersic = m.add_profile("sersic");
		sersic->parameter("xcen", 50.0);
		sersic->parameter("ycen", 50.0);
		sersic->parameter("mag", 13.0);
		sersic->parameter("re", 5.0);
		sersic->parameter("nser", 8.0);
		sersic->parameter("ang", 80.);
		sersic->parameter("axrat", 0.6);
		sersic->parameter("box", -0.3);
		_check_images_within_tolerance(m);
	}

	void test_convolver() {
		if( !openCLFixtures.opencl_env ) {
			TS_SKIP("No OpenCL environment found to run OpenCL tests");
		}
		_check_convolver(OpenCLConvolver{openCLFixtures.opencl_env});

		}
	}

};

#endif // PROFIT_OPENCL