/**
 * library tests
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2018
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

#define PROFIT_TEST_NO_LIBRARY_FIXTURE
#include "common_test_setup.h"

class TestLibrary : public CxxTest::TestSuite {

public:
	void test_init_finish() {
		// simply run them and make sure they don't fail
		TSM_ASSERT("Failed to initialize library", profit::init());
		TS_ASSERT(profit::init_diagnose().empty())
		profit::finish();
		TS_ASSERT(profit::finish_diagnose().empty())
	}

};