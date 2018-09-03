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

#include <vector>

#include "common_test_setup.h"

using namespace profit;

std::vector<const char *> all_radial = {
	"brokenexp",
	"coresersic",
	"ferrer",
	"king",
	"moffat",
	"sersic"
};

class TestRadial : public CxxTest::TestSuite {

public:

	void test_create_default(void) {
		for(auto pname: all_radial) {

			// 100x100 model, with profile centered at (50,50)
			Model m {100, 100};
			auto radialp = m.add_profile(pname);
			radialp->parameter("xcen", 50.);
			radialp->parameter("ycen", 50.);

			// we don't assert anything yet, only check that the profile can be
			// constructed successfully with its default values
			m.evaluate();
		}
	}

	void test_create_boxy(void) {
		for(auto pname: all_radial) {

			// 100x100 model, with profile centered at (50,50), box=0.5
			Model m {100, 100};
			auto radialp = m.add_profile(pname);
			radialp->parameter("xcen", 50.);
			radialp->parameter("ycen", 50.);
			radialp->parameter("box", 0.1);

			// we don't assert anything yet, only check that the profile can be
			// constructed successfully with its default values
			m.evaluate();
		}
	}

	void test_calcmask(void) {

		Model m {3, 3};
		m.set_psf({{1,1,1,1}, 2, 2});
		m.set_magzero(0);

		auto sp = m.add_profile("sersic");
		sp->parameter("xcen", 1.);
		sp->parameter("ycen", 1.);
		sp->parameter("re", 10.);
		sp->parameter("rscale_max", 10.);
		sp->parameter("mag", 0.);
		sp->parameter("adjust", false);

		/* Some on, some off */
		Mask mask {{false, true, true,
		            false, true, false,
		            true, true, false}, 3, 3};
		m.set_mask(mask);
		auto image = m.evaluate();
		auto image_it = image.begin();
		auto mask_it = mask.begin();
		for(; image_it != image.end(); image_it++, mask_it++) {
			if( *mask_it ) {
				TS_ASSERT_DIFFERS(0, *image_it);
			}
			else {
				TS_ASSERT_EQUALS(0, *image_it);
			}
		}

		/* All on */
		mask = {{true, true, true,
		         true, true, true,
		         true, true, true}, 3, 3};
		m.set_mask(mask);
		image = m.evaluate();
		for(auto pixel: image) {
			TS_ASSERT_DIFFERS(0, pixel);
		}

		/* No mask */
		m.set_mask(Mask {});
		image = m.evaluate();
		for(auto pixel: image) {
			TS_ASSERT_DIFFERS(0, pixel);
		}
	}

};
