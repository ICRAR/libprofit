/**
 * utils tests
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


#include <cmath>
#include <cstdlib>
#include <iostream>

#include <cxxtest/TestSuite.h>

#include "profit/config.h"
#include "profit/utils.h"

#ifdef PROFIT_USES_GSL
	#include <gsl/gsl_errno.h>
#endif


using namespace profit;

class TestUtils : public CxxTest::TestSuite {

public:

#ifdef PROFIT_USES_GSL
	void setUp() {
		gsl_set_error_handler_off();
	}
#endif

	void test_gamma(void) {
		// not defined for negative integers and 0
		for(double x: {-100., -50., -1., 0.}) {
			TS_ASSERT(std::isnan(gammafn(x)));
		}
		// normally approximates to 0 for large negative numbers (non-integer)
		for(double x: {-1000.1, -2000.5}) {
			double result = gammafn(x);
			TS_ASSERT(result == 0.);
		}
		// defined for values up to ~171
		for(double x: {-100.1, -50.1, -1.1, 0.1, 1., 2., 3.1, 10., 50., 70., 100.}) {
			double result = gammafn(x);
			TS_ASSERT(!std::isnan(result) && !std::isinf(result));
		}
		// and returns +Inf for big positives
		for(double x: {200., 2500., 10000.}) {
			TS_ASSERT(std::isinf(gammafn(x)));
		}
	}

	void test_beta() {

		// not defined for a or b < 0
		for(double x: {-100., -10., -5., -4.5, -1., -0.1}) {
			double result1 = beta(1, x);
			double result2 = beta(x, 1);
			TS_ASSERT(std::isnan(result1));
			TS_ASSERT(std::isnan(result2));
		}

		// returns +Inf for a or b == 0
		TS_ASSERT(std::isinf(beta(0, 1)));
		TS_ASSERT(std::isinf(beta(1, 0)));

		// normal values...
		for(double a: {1., 2., 3., 4., 10., 100., 200.}) {
			for(double b: {1., 2., 3., 4., 10., 100., 200.}) {
				double result = beta(a,b);
				TS_ASSERT(!std::isinf(result));
				TS_ASSERT(!std::isnan(result));
				TS_ASSERT_DIFFERS(0., result)
			}
		}

		// for big values 0 is returned
		for(double x: {1000., 2000., 3000., 3000.}) {
			TS_ASSERT_EQUALS(0., beta(x,x));
		}
	}

	void test_get_profit_home() {

		// Test both the normal, HOME-based profit home directory
		// and the environment variable-based directory
		auto run_test = []() {
			auto home = get_profit_home();
			TSM_ASSERT_DIFFERS("Profit's home is not empty", "", home);
			TS_ASSERT(dir_exists(home));
		};

		// Small utility to keep the previous PROFIT_HOME, if any
		struct env_keeper {
			std::string _home;
			env_keeper() : _home() {
				auto h = std::getenv("PROFIT_HOME");
				if (h != nullptr) {
					_home = h;
				}
			}
			~env_keeper() {
				profit::setenv("PROFIT_HOME", _home);
			}
		} keeper;

		// Run the tests both with the default profit home and with a hardcoded one
		run_test();
		profit::setenv("PROFIT_HOME", ".profit");
		run_test();
		recursive_remove(".profit");
	}

	void test_recursive_remove() {
		auto home = get_profit_home();
		recursive_remove(home);
		TS_ASSERT(!dir_exists(home));
	}

	void test_recursive_remove_failures() {
		TS_ASSERT_THROWS(recursive_remove(".this_shouldn_exist_at_all"), std::runtime_error);
	}

};