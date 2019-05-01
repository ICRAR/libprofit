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

#include <random>
#include "common_test_setup.h"

using namespace profit;

class TestConvolver : public CxxTest::TestSuite {

private:

	Image uniform_random_image(Dimensions dim)
	{
		Image image(dim);
		std::random_device dev;
		std::default_random_engine engine(dev());
		std::uniform_real_distribution<double> uniform(0, 1);
		for(auto &d: image) {
			d = uniform(engine);
		}
		return image;
	}

	void _pixels_within_tolerance(const Image &original_im, const Image &new_im,
	                              unsigned int i, double tolerance) {

		unsigned int width = original_im.getWidth();
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

	void _check_convolver(ConvolverPtr &&otherConvolver) {
		for(auto im_dim: {100U, 101U}) {
			for(auto krn_dim: {24U, 25U}) {
				Mask mask;
				auto src = uniform_random_image({im_dim, im_dim});
				auto krn = uniform_random_image({krn_dim, krn_dim});
				auto bConvolver = create_convolver(ConvolverType::BRUTE_OLD);
				Image result1 = bConvolver->convolve(src, krn, mask);
				Image result2 = otherConvolver->convolve(src, krn, mask);
				for(unsigned int i = 0; i < src.size(); i++) {
					// Hopefully within 0.1% of error?
					_pixels_within_tolerance(result1, result2, i, 1e-3);
				}
			}
		}
	}

	void _test_openmp_convolver(ConvolverType type) {

		auto src = uniform_random_image({50, 50});
		auto krn = uniform_random_image({25, 25});
		Mask mask;

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

	void _test_masked_convolution(ConvolverType type)
	{
		// The mask kind of draws an arrow pointing to the upper-left corner
		std::vector<double> src {1, 1, 1, 1,
		                         1, 1, 1, 1,
		                         1, 1, 1, 1,
		                         1, 1, 1, 1};
		std::vector<bool> mask {true, true, true, false,
		                        true, true, false, false,
		                        true, false, true, false,
		                        false, false, false, true};
		std::vector<double> krn {1, 1, 1,
		                         1, 1, 1,
		                         1, 1, 1};

		Image src_img(src, 4, 4);
		Image krn_img(krn, 3, 3);
		Mask m(mask, 4, 4);

		auto convolver = create_convolver(type, ConvolverCreationPreferences());
		auto result = convolver->convolve(src_img, krn_img, m);

		// Only the pixels where the mask is true should be set; the rest should
		// be zero
		auto res_it = result.begin();
		auto mask_it = m.begin();
		size_t i = 0;
		for(; res_it != result.end(); res_it++, mask_it++, i++) {
			if (*mask_it) {
				std::ostringstream msg;
				msg << "Pixel [" << (i % 4) << "," << (i / 4) << "] is zero, but should not be";
				TSM_ASSERT_DIFFERS(msg.str(), 0, *res_it);
			}
			else {
				std::ostringstream msg;
				msg << "Pixel [" << (i % 4) << "," << (i / 4) << "] is not zero, but should";
				TSM_ASSERT_EQUALS(msg.str(), 0, *res_it);
			}
		}
	}

	void _test_simd_convolver(simd_instruction_set instruction_set)
	{
		std::ostringstream msg;
		msg << instruction_set << " should be invalid";
		ConvolverCreationPreferences prefs;
		prefs.instruction_set = instruction_set;
		if (has_simd_instruction_set(instruction_set)) {
			_check_convolver(create_convolver(ConvolverType::BRUTE, prefs));
		}
		else {
			TSM_ASSERT_THROWS(msg.str(), create_convolver(ConvolverType::BRUTE, prefs), const invalid_parameter &);
		}
	}

public:

	void test_new_bruteforce_convolver() {
		ConvolverCreationPreferences prefs;
		_check_convolver(create_convolver(ConvolverType::BRUTE, prefs));
	}

	void test_brute_convolver_AVX() {
		_test_simd_convolver(simd_instruction_set::AVX);
	}

	void test_brute_convolver_SSE2() {
		_test_simd_convolver(simd_instruction_set::SSE2);
	}

	void test_brute_convolver_NONE() {
		_test_simd_convolver(simd_instruction_set::NONE);
	}

	void test_brute_convolver_AUTO() {
		_test_simd_convolver(simd_instruction_set::AUTO);
	}

	void test_old_bruteforce_openmp() {
		_test_openmp_convolver(ConvolverType::BRUTE_OLD);
	}

	void test_new_bruteforce_openmp() {
		_test_openmp_convolver(ConvolverType::BRUTE);
	}

	void test_masked_convolution() {
		_test_masked_convolution(ConvolverType::BRUTE_OLD);
		_test_masked_convolution(ConvolverType::BRUTE);
	}
};