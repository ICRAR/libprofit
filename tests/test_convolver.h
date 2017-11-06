//
// Convolver-specific tests
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
// Copyright by UWA (in the framework of the ICRAR)
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307  USA
//

#include <cxxtest/TestSuite.h>

#include "profit/profit.h"

using namespace profit;

class TestConvolver : public CxxTest::TestSuite {

private:

	void _pixels_within_tolerance(std::vector<double> original_im, std::vector<double> new_im,
	                              unsigned int i, unsigned int width,
	                              double tolerance) {

		double original = original_im[i];
		double new_pixel = new_im[i];
		auto diff = std::abs(original - new_pixel);
		if ( !diff ) {
			// all good
			return;
		}

		// avoid NaNs due to divide-by-zero
		auto denomin = original;
		if ( !denomin ) {
			denomin = new_pixel;
		}
		auto relative_diff = diff / denomin;

		std::ostringstream msg;
		msg << "Pixel [" << i % width << "," << i / width << "] has values that are too different: ";
		msg << original << " v/s " << new_pixel;
		TSM_ASSERT_LESS_THAN_EQUALS(msg.str(), relative_diff, tolerance);
	}

	void _check_convolver(Convolver &&otherConvolver) {
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

				auto bConvolver = create_convolver(ConvolverType::BRUTE_OLD);
				Image result1 = bConvolver->convolve(src, krn, mask);
				Image result2 = otherConvolver.convolve(src, krn, mask);
				for(unsigned int i = 0; i < src.getSize(); i++) {
					// Hopefully within 0.1% of error?
					_pixels_within_tolerance(result1.getData(), result2.getData(), i, src.getWidth(), 1e-12);
				}
			}
		}
	}

	void _test_openmp_convolver(ConvolverType type) {

		// Random images
		auto im_dim = 50;
		auto krn_dim = 25;
		Image src(im_dim, im_dim);
		Image krn(krn_dim, krn_dim);
		Mask mask;
		for(auto &d: src.getData()) {
			d = (rand() % 10000) / 10000.0;
		}
		for(auto &d: krn.getData()) {
			d = (rand() % 10000) / 10000.0;
		}

		// A normal and an OpenMP-accelerated brute-force convolver
		auto prefs = ConvolverCreationPreferences();
		prefs.omp_threads = 2;
		auto brute = create_convolver(type);
		auto openmp = create_convolver(type, prefs);

		// Results should be fully identical
		auto res1 = brute->convolve(src, krn, mask);
		auto res2 = openmp->convolve(src, krn, mask);
		TS_ASSERT(res1 == res2);

	}

public:

	void test_new_bruteforce_convolver() {
		_check_convolver(AssociativeBruteForceConvolver(1));
	}

	void test_old_bruteforce_openmp() {
		_test_openmp_convolver(ConvolverType::BRUTE_OLD);
	}

	void test_new_bruteforce_openmp() {
		_test_openmp_convolver(ConvolverType::BRUTE);
	}

};