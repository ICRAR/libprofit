/**
 * FFT convolution tests
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

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <memory>
#include <string>
#include <vector>

#include <cxxtest/GlobalFixture.h>
#include <cxxtest/TestSuite.h>

#include "profit/profit.h"

using namespace profit;

/**
 * Create and keep a pointer to an appropriate plan for the given sizes
 */
class FFTFixtures : CxxTest::GlobalFixture {

public:
	FFTFixtures(const Dimensions &dims, const Dimensions &psf_dims, bool reuse_psf_fft) :
		tolerance(0.01),
		fft_convolver(nullptr),
		psf(psf_dims),
		dims(dims)
	{

		if (not has_fftw()) {
			return;
		}

		ConvolverCreationPreferences prefs;
		prefs.src_dims = dims;
		prefs.krn_dims = psf_dims;
		prefs.effort = effort_t::ESTIMATE;
		prefs.omp_threads = 1;
		prefs.reuse_krn_fft = reuse_psf_fft;
		fft_convolver = create_convolver(ConvolverType::FFT, prefs);
		// a random psf
		unsigned int seed = (unsigned int)time(NULL);
		rand_r(&seed);
		for (unsigned int i = 0; i < psf.size(); i++) {
			psf[i] = (rand() % 10000) / 10000.0;
		}

		// For manual checks, if wanted
		if( const char *tol = std::getenv("LIBPROFIT_FFT_TOLERANCE") ) {
			tolerance = std::stod(tol);
		}
	}

	double tolerance;
	ConvolverPtr fft_convolver;
	Image psf;
	Dimensions dims;
};

// all combinations of: even/odd image, even/odd psf, psf fft reuse/no_reuse
static FFTFixtures fixtures[] = {
	{{100, 100}, {25, 25}, true},
	{{100, 100}, {26, 26}, true},
	{{99, 99}, {25, 25}, true},
	{{99, 99}, {26, 26}, true},
	{{100, 100}, {25, 25}, false},
	{{100, 100}, {26, 26}, false},
	{{99, 99}, {25, 25}, false},
	{{99, 99}, {26, 26}, false}
};

class TestFFT : public CxxTest::TestSuite {

public:

	void _check_fftw_support() {
		if (not has_fftw()) {
			TS_SKIP("No FFTW support found, cannot run this test");
		}
	}

	void _check_images_within_tolerance(Model &m) {
		for(auto &fixture: fixtures) {
			_check_images_within_tolerance_with_fixtures(m, fixture);
		}
	}

	void _check_images_within_tolerance_with_fixtures(Model &m, FFTFixtures &fftFixtures)
	{

		if( !fftFixtures.fft_convolver ) {
			TS_SKIP("No FFTPlan found to run FFT tests with this fixture");
		}

		auto width = fftFixtures.dims.x;
		m.set_dimensions(fftFixtures.dims);
		m.set_psf(fftFixtures.psf);

		// evaluate normally first, and then using the FFTPlan
		m.set_convolver(nullptr);
		std::vector<double> original = m.evaluate();
		m.set_convolver(fftFixtures.fft_convolver);
		std::vector<double> fft_produced = m.evaluate();

		// Pixel by pixel the images should be fairly similar
		for(unsigned int i=0; i!=original.size(); i++) {

			auto original_pixel = original[i];
			auto fft_pixel = fft_produced[i];

			auto diff = std::abs(original_pixel - fft_pixel);
			if ( !diff ) {
				// all good
				continue;
			}

			// avoid NaNs due to divide-by-zero
			// Also, when we have zero on a pixel with simple brute-force convolution
			// it will be very hard to have a zero on the FFT-based convolution.
			// The only thing we can really do is to assert that the value is
			// indeed very low compared to the rest of the image
			auto denomin = original_pixel;
			if ( !denomin ) {
				continue;
			}

			std::ostringstream msg;
			auto relative_diff = diff / denomin;
			msg << "Pixel [" << i % width << "," << i / width << "] has values that are too different: ";
			msg << original_pixel << " v/s " << fft_pixel;
			TSM_ASSERT_LESS_THAN_EQUALS(msg.str(), relative_diff, fftFixtures.tolerance);
		}
	}

	void test_fft_brokenexp() {
		_check_fftw_support();
		Model m;
		auto brokenexp = m.add_profile("brokenexp");
		brokenexp->parameter("convolve", true);
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

	void test_fft_coresersic() {
		_check_fftw_support();
		Model m;
		auto coresersic = m.add_profile("coresersic");
		coresersic->parameter("convolve", true);
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

	void test_fft_ferrer() {
		_check_fftw_support();
		Model m;
		auto ferrer = m.add_profile("ferrer");
		ferrer->parameter("convolve", true);
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

	void test_fft_king() {
		_check_fftw_support();
		Model m;
		auto king = m.add_profile("king");
		king->parameter("convolve", true);
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

	void test_fft_moffat() {
		_check_fftw_support();
		Model m;
		auto moffat = m.add_profile("moffat");
		moffat->parameter("convolve", true);
		moffat->parameter("xcen", 34.);
		moffat->parameter("ycen", 74.);
		moffat->parameter("mag", 10.);
		moffat->parameter("fwhm", 10.0);
		moffat->parameter("con", 10.);
		moffat->parameter("ang", 45.);
		moffat->parameter("axrat", 0.95);
		moffat->parameter("box", 0.1);
		_check_images_within_tolerance(m);
	}

	void test_fft_sersic() {
		_check_fftw_support();
		Model m;
		auto sersic = m.add_profile("sersic");
		sersic->parameter("convolve", true);
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

	void test_fftconvolver_reusage() {
		_check_fftw_support();
		Mask mask;
		Image src(100, 100);
		Image krn(25, 25);
		for(auto &d: src) {
			d = (rand() % 10000) / 10000.0;
		}
		for(auto &d: krn) {
			d = (rand() % 10000) / 10000.0;
		}

		ConvolverCreationPreferences prefs;
		prefs.src_dims = {100, 100};
		prefs.krn_dims = {25, 25};
		prefs.effort = effort_t::ESTIMATE;
		prefs.reuse_krn_fft = true;
		auto convolver = create_convolver(ConvolverType::FFT, prefs);
		Image result1 = convolver->convolve(src, krn, mask, false);
		Image result2 = convolver->convolve(src, krn, mask, false);
		Image result3 = convolver->convolve(src, krn, mask, false);
		TS_ASSERT(result1 == result2);
		TS_ASSERT(result2 == result3);
	}

};