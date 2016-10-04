/**
 * radial profile tests
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

#include "profit/sersic.h"

using namespace profit;

class TestRadial : public CxxTest::TestSuite {

public:

	void test_calcmask(void) {

		Model m;
		m.width = 3;
		m.height = 3;
		m.psf_width = 2;
		m.psf_height = 2;
		m.magzero = 0;
		m.psf = {1,1,1,1};

		auto &psfprof = static_cast<SersicProfile &>(m.add_profile("sersic"));
		psfprof.xcen = 1;
		psfprof.ycen = 1;
		psfprof.re = 10;
		psfprof.rscale_max = 10;
		psfprof.mag = 0;
		psfprof.adjust = false;

		/* Some on, some off */
		m.calcmask = {false, true, true,
		              false, true, false,
		              true, true, false};
		std::vector<double> image = m.evaluate();
		for(auto j=0; j!=m.width; j++) {
			for(auto i=0; i!=m.height; i++) {
				auto idx = i + j*m.width;
				if( m.calcmask[idx] ) {
					TS_ASSERT_DIFFERS(0, image[idx]);
				}
				else {
					TS_ASSERT_EQUALS(0, image[idx]);
				}
			}
		}

		/* All on */
		m.calcmask = {true, true, true,
		              true, true, true,
		              true, true, true};
		image = m.evaluate();
		for(auto j=0; j!=m.width; j++) {
			for(auto i=0; i!=m.height; i++) {
				auto idx = i + j*m.width;
				TS_ASSERT_DIFFERS(0, image[idx]);
			}
		}

		/* No mask */
		m.calcmask = {};
		image = m.evaluate();
		for(auto j=0; j!=m.width; j++) {
			for(auto i=0; i!=m.height; i++) {
				auto idx = i + j*m.width;
				TS_ASSERT_DIFFERS(0, image[idx]);
			}
		}
	}

};
