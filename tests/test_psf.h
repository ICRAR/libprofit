/**
 * psf profile tests
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


#include "common_test_setup.h"

using namespace profit;

class TestPsf : public CxxTest::TestSuite {

public:

	void test_exact_pixels(void) {

		Model m {10, 10};
		m.set_psf({{1,1,1,1}, 2, 2});
		m.set_magzero(0);

		auto psfprof = m.add_profile("psf");
		psfprof->parameter("xcen", 2.);
		psfprof->parameter("ycen", 2.);
		psfprof->parameter("mag", 0.);

		auto image = m.evaluate();
		for(auto &pixel: image) {
			TS_ASSERT_DELTA(0.25, pixel, 1e9);
		}

	}

};