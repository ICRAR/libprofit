/**
 * profile-wide tests
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

#include <utility>

#include "common_test_setup.h"

using namespace profit;

class TestProfile : public CxxTest::TestSuite {

private:

	void test_params_i(const std::string &profile_name,
	                   const std::vector<std::string> &unknown_names,
	                   const std::vector<std::string> &bool_names,
	                   const std::vector<std::string> &uint_names,
	                   const std::vector<std::string> &double_names) {

		profit::Model m;
		auto prof = m.add_profile(profile_name);

		// invalid parameters
		for(auto name: unknown_names) {
			TS_ASSERT_THROWS(prof->parameter(name, true), const invalid_parameter &);
			TS_ASSERT_THROWS(prof->parameter(name, 1u), const invalid_parameter &);
			TS_ASSERT_THROWS(prof->parameter(name, 1.), const invalid_parameter &);
		}

		// valid boolean parameters
		for(auto name: bool_names) {
			prof->parameter(name, true);
			TS_ASSERT_THROWS(prof->parameter(name, 1u), const invalid_parameter &);
			TS_ASSERT_THROWS(prof->parameter(name, 1.), const invalid_parameter &);
		}

		// valid uint parameters
		for(auto name: uint_names) {
			TS_ASSERT_THROWS(prof->parameter(name, true), const invalid_parameter &);
			prof->parameter(name, 1u);
			TS_ASSERT_THROWS(prof->parameter(name, 1.), const invalid_parameter &);
		}

		// valid double parameters
		for(auto name: double_names) {
			TS_ASSERT_THROWS(prof->parameter(name, true), const invalid_parameter &);
			TS_ASSERT_THROWS(prof->parameter(name, 1u), const invalid_parameter &);
			prof->parameter(name, 1.);
		}

	}

	void test_param_values(const std::string &profile_name,
	                       const std::string &param_name,
	                       const std::vector<double> &allowed_values,
	                       const std::vector<double> &invalid_values,
	                       const std::vector<std::pair<std::string, double>> &fixed_vals = std::vector<std::pair<std::string, double>>()) {

		Model m {10, 10};
		m.set_dry_run(true);
		auto p = m.add_profile(profile_name);

		// check that by default all profiles generate valid images
		m.evaluate();

		// set the fixed vals before any any param_name test
		for(auto &fixed_val: fixed_vals) {
			p->parameter(std::get<0>(fixed_val), std::get<1>(fixed_val));
		}

		// allowed values don't throw exceptions
		for(auto v: allowed_values) {
			p->parameter(param_name, v);
			m.evaluate();
		}

		// invalid values do
		for(auto v: invalid_values) {
			p->parameter(param_name, v);
			TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		}
	}

	void test_param_positive(const std::string &profile_name,
	                         const std::string &param_name) {
		test_param_values(profile_name, param_name, {0.1, 0.2, 1, 4}, {-20, -10, -1, -0.001, 0});
	}

	void test_param_positive_and_zero(const std::string &profile_name,
	                                  const std::string &param_name) {
		test_param_values(profile_name, param_name, {0, 0.1, 0.2, 1, 4}, {-20, -10, -1, -0.001});
	}

	void test_radial_params(const std::string &profile_name) {
		test_params_i(profile_name, {"Xcen", "yCen", "magnitude", "axisrat"},
		                            {"adjust", "rough"}, {"max_recursions", "resolution"},
		                            {"xcen", "ycen", "mag", "ang", "axrat", "box", "acc", "rscale_switch", "rscale_max"});
		test_param_values(profile_name, "axrat",
		                  {0.1, 0.2, 0.5, 0.99, 1},
		                  {-1000, -100, -10, -1, -0.1, 0, 1.001, 2, 4, 8, 16, 1000});
		test_param_values(profile_name, "box",
		                  {-1.5, -1, 0, 0.5},
		                  {-100, -10, -5, -2.001, -2});
	}

public:

	void test_brokenexp_parameters() {
		test_radial_params("brokenexp");
		test_params_i("brokenexp",
				{"b", "H1", "H2", "h3", "1h", "unknown"},
				{},
				{},
				{"a", "h1", "h2", "rb"});
		test_param_positive("brokenexp", "rb");
		// h2 should be <= h1, so whenever we test h2 we should also set h1 and viceversa
		test_param_values("brokenexp", "h1", {0.1, 1, 5, 10}, {-10, -5, -2, -1, 0}, {std::make_pair("h2", 0.09)});
		test_param_values("brokenexp", "h2", {0.1, 1, 5, 10}, {-10, -5, -2, -1, 0}, {std::make_pair("h1", 1000.)});
	}

	void test_coresersic_parameters() {
		test_radial_params("coresersic");
		test_params_i("coresersic",
				{"v", "c", "Nser", "Re", "RE", "unknown"},
				{},
				{},
				{"a", "b", "nser", "rb", "re"});

		test_param_positive("coresersic", "re");
		test_param_positive("coresersic", "rb");
		test_param_positive("coresersic", "nser");
		test_param_positive("coresersic", "a");
		test_param_values("coresersic", "b", {-2, -1, 0, 1, 1.8, 1.9}, {2, 3, 4, 10});
	}

	void test_ferrer_parameters() {
		test_radial_params("ferrer");
		test_params_i("ferrer",
				{"v", "c", "rin", "Rout", "unknown"},
				{},
				{},
				{"a", "b", "rout"});
		test_param_positive("ferrer", "rout");
		test_param_positive_and_zero("ferrer", "a");
		test_param_values("ferrer", "b", {-2, -1, 0, 1, 1.8, 1.9, 2}, {2.1, 3, 4, 10});
	}

	void test_king_parameters() {
		test_radial_params("king");
		test_params_i("king",
				{"b", "c", "r", "rg", "unknown"},
				{},
				{},
				{"a", "rc", "rt"});
		test_param_positive("king", "rc");
		test_param_positive("king", "rt");
		test_param_positive_and_zero("king", "a");
	}

	void test_moffat_parameters() {
		test_radial_params("moffat");
		test_params_i("moffat",
				{"fmwh", "fwmh", "fmhw", "unknown"},
				{},
				{},
				{"con", "fwhm"});
		test_param_positive("moffat", "fwhm");
		test_param_positive_and_zero("moffat", "con");
	}

	void test_psf_parameters() {
		test_params_i("psf",
				{"mag1", "YCEN", "xcenter"},
				{},
				{},
				{"xcen", "ycen", "mag"});
	}

	void test_sky_parameters() {
		test_params_i("sky",
				{"bg2", "bg1", "BG", "unknown"},
				{},
				{},
				{"bg"});
	}

	void test_sersic_parameters() {
		test_radial_params("sersic");
		test_params_i("sersic",
				{"unknown", "RE", "nser "},
				{"rescale_flux"},
				{},
				{"nser", "re"});
		test_param_positive("sersic", "nser");
		test_param_positive("sersic", "re");
	}

};