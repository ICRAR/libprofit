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
#include "profit/sky.h"

using namespace profit;

class TestPsf : public CxxTest::TestSuite {

public:

	void test_valid_dimensions(void) {

		Model m;

		// Only the final combination is valid
		m.width = 0; m.height = 0;
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		m.width = 0; m.height = 1;
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		m.width = 1; m.height = 0;
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		m.width = 1; m.height = 1;
		m.evaluate(); // fine...

		// We generate the correct size
		m.width = 100;
		m.height = 100;
		TS_ASSERT_EQUALS(m.width * m.height, m.evaluate().size())
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
		TS_ASSERT_EQUALS(false, m.has_profiles());

		// these are fine
		for(auto name: {"brokenexp", "coresersic", "ferrer", "ferrers",
		                "king", "moffat", "psf", "sersic", "sky"}) {
			m.add_profile(name);
		}
		TS_ASSERT_EQUALS(true, m.has_profiles());

	}

	void test_valid_psf(void) {

		Model m;
		m.height = m.width = 2;

		// so far so good...
		m.evaluate();

		// add a profile, ask to convolve it, but don't give a psf
		auto &skyp = static_cast<SkyProfile &>(m.add_profile("sky"));
		skyp.convolve = true;
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);

		// Give a psf, but still without dimensions information
		m.psf = {1,1};
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);

		// 0 dimensions
		m.psf_height = 0; m.psf_width = 0;
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		m.psf_height = 0; m.psf_width = 1;
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		m.psf_height = 1; m.psf_width = 0;
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		// This is still invalid because the PSF vector has 2 elements
		m.psf_height = 1; m.psf_width = 1;
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);

		// finally everything is in order...
		m.psf_height = 1; m.psf_height = 2;
		m.evaluate();
	}
};