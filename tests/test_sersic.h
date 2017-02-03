/**
 * sersic profile tests
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

#include <vector>

#include <cxxtest/TestSuite.h>

#include "profit/profit.h"

using namespace profit;

class TestSersic : public CxxTest::TestSuite {

private:

	void _test_different_nser(double box) {

		// 50x50 model, with profile centered at (25,25)
		Model m;
		m.width = m.height = 50;
		auto sersicp = m.add_profile("sersic");
		sersicp->parameter("xcen", 25.);
		sersicp->parameter("ycen", 25.);
		sersicp->parameter("box", box);

		for(auto nser: {0.5, 1., 1.5, 2., 3., 4., 8., 16.}) {
			// we don't assert anything yet, only check that the profile can be
			// constructed successfully with the given values
			sersicp->parameter("nser", nser);
			m.evaluate();
		}

	}

public:

	void test_different_nser(void) {
		_test_different_nser(0);
	}

	void test_different_nser_boxy(void) {
		_test_different_nser(0.1);
	}

};
