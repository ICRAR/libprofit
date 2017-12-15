/**
 * sky profile tests
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

class TestSky : public CxxTest::TestSuite {

public:

	void test_bg(void) {

		Model m {2, 2};

		auto skyprof = m.add_profile("sky");
		skyprof->parameter("bg", 1.);

		std::vector<double> image = m.evaluate();
		for(auto pixel: image) {
			TS_ASSERT_EQUALS(1, pixel);
		}

	}

	void test_with_mask(void) {

		Model m {2, 2};
		m.set_mask({{true, true, true, false}, 2, 2});

		auto skyprof = m.add_profile("sky");
		skyprof->parameter("bg", 5.);

		auto image = m.evaluate();
		for(int idx: {0,1,2}) {
			TS_ASSERT_EQUALS(5, image[idx]);
		}
		TS_ASSERT_EQUALS(0, image[3]);

	}
};