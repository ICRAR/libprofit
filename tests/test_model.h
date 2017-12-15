/**
 * model tests
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

#include <cxxtest/TestSuite.h>

#include "profit/profit.h"

using namespace profit;

class TestModel : public CxxTest::TestSuite {

public:

	void test_valid_dimensions(void) {

		Model m;

		// Only the final combination is valid
		m.set_dimensions({0, 0});
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		m.set_dimensions({0, 1});
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		m.set_dimensions({1, 0});
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		m.set_dimensions({1, 1});
		m.evaluate(); // fine...

		// We generate the correct size
		Dimensions dims {100, 100};
		m.set_dimensions(dims);
		TS_ASSERT_EQUALS(dims.x * dims.y, m.evaluate().size());
	}

	void test_valid_scales(void) {

		Model m {1, 1};

		for(auto scale: {-2, -1, 0}){
			m.set_image_pixel_scale({scale, 1});
			TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		}
		m.set_image_pixel_scale({0.1, 1});
		m.evaluate(); // fine

		for(auto scale: {-2, -1, 0}){
			m.set_image_pixel_scale({0.1, scale});
			TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		}
		m.set_image_pixel_scale({0.1, 0.1});
		m.evaluate(); // fine

	}

	void test_valid_profiles(void) {

		Model m;
		for(auto name: {"unknown", "sersic1", "Sersic", "sersi", " sersic", "sersic "}) {
			TS_ASSERT_THROWS(m.add_profile(name), const invalid_parameter &);
		}
		TS_ASSERT_EQUALS(false, m.has_profiles());

		// these are fine
		for(auto name: {"brokenexp", "coresersic", "ferrer", "ferrers",
		                "king", "moffat", "psf", "sersic", "sky"}) {
			m.add_profile(name);
		}
		TS_ASSERT_EQUALS(true, m.has_profiles());

	}

	void test_valid_psf(void) {

		Model m {2, 2};

		// so far so good...
		m.evaluate();

		// add a profile, ask to convolve it, but don't give a psf
		auto skyp = m.add_profile("sky");
		skyp->parameter("convolve", true);
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);

		// Give a psf, everything is in order...
		m.set_psf({{1, 1}, 1, 2});
		m.evaluate();
	}

	void _add_sersic(Model &m, double xcen, double ycen, double re, bool convolve = false) {
		auto sersic = m.add_profile("sersic");
		sersic->parameter("xcen", xcen);
		sersic->parameter("ycen", ycen);
		sersic->parameter("re", re);
		sersic->parameter("convolve", convolve);
	}

	void test_profile_images_addition() {

		// three individual model images are summed up
		Model m1(100, 100);
		_add_sersic(m1, 50, 50, 10);
		auto image1 = m1.evaluate();

		Model m2(100, 100);
		_add_sersic(m2, 30, 10, 16);
		auto image2 = m2.evaluate();

		Model m3(100, 100);
		_add_sersic(m3, 23, 89, 1.2);
		auto image3 = m3.evaluate();

		// image1 holds the final result
		image1 = image1 + image2 + image3;

		// A single model image with all profile images
		Model m4(100, 100);
		_add_sersic(m4, 50, 50, 10);
		_add_sersic(m4, 30, 10, 16);
		_add_sersic(m4, 23, 89, 1.2);
		auto image4 = m4.evaluate();

		// They should be the same! We add them in the same order to make sure
		// that floating-point rounding yields the same result
		TS_ASSERT(image1 == image4);
	}

	void test_profile_images_addition_after_convolving() {

		auto convolver = create_convolver(ConvolverType::BRUTE);
		auto psf = {0., 1., 2., 3.};

		// two individual model images are summed up
		// the second is actually convolved with a psf
		Model m1(100, 100);
		_add_sersic(m1, 50, 50, 10);
		auto image1 = m1.evaluate();

		Model m2(100, 100);
		m2.set_convolver(convolver);
		m2.set_psf({psf, 2, 2});
		_add_sersic(m2, 30, 10, 16, true);
		auto image2 = m2.evaluate();

		// image1 holds the final result
		image1 = image1 + image2;

		// A single model image with all profile images
		Model m3(100, 100);
		m3.set_convolver(convolver);
		m3.set_psf({psf, 2, 2});
		_add_sersic(m3, 50, 50, 10);
		_add_sersic(m3, 30, 10, 16, true);
		auto image3 = m3.evaluate();

		// They should be the same! We add them in the same order to make sure
		// that floating-point rounding yields the same result
		TS_ASSERT(image1 == image3);

	}
};