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

class TestPsf : public CxxTest::TestSuite {

public:

	void test_valid_dimensions(void) {

		Model m;

		m.width = 0;
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		m.height = 0;
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);

		m.width = 1;
		m.height = 1;
		m.evaluate(); // fine...

	}

	void test_valid_scales(void) {

		Model m;
		m.width = 1;
		m.height = 1;

		m.scale_y = 1;
		for(auto scale: {-2, -1, 0}){
			m.scale_x = scale;
			TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		}
		m.scale_x = 0.1;
		m.evaluate(); // fine

		for(auto scale: {-2, -1, 0}){
			m.scale_y = scale;
			TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		}
		m.scale_y = 0.1;
		m.evaluate(); // fine

	}

	void test_valid_profiles(void) {

		Model m;
		for(auto name: {"unknown", "sersic1", "Sersic", "sersi", " sersic", "sersic "}) {
			TS_ASSERT_THROWS(m.add_profile(name), const invalid_parameter &);
		}

		// these are fine
		for(auto name: {"brokenexp", "coresersic", "ferrer", "king", "moffat", "psf", "sersic", "sky"}) {
			m.add_profile(name);
		}

	}

};