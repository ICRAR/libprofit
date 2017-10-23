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

public:

	void test_openmp_convolver() {

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
		auto brute = create_convolver(ConvolverType::BRUTE);
		auto openmp = create_convolver(ConvolverType::BRUTE, prefs);

		// Results should be fully identical
		auto res1 = brute->convolve(src, krn, mask);
		auto res2 = openmp->convolve(src, krn, mask);
		TS_ASSERT(res1 == res2);

	}

};