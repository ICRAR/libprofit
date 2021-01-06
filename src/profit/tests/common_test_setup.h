//
// Common test functionality
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

#include <cxxtest/GlobalFixture.h>
#include <cxxtest/TestSuite.h>

#include "profit/profit.h"


// Nice printing for some libprofit objects
namespace CxxTest {

// base class to be used with all classes supporting operator<<
template <typename T>
class StreamedString
{
	std::string as_string;
public:
	explicit StreamedString(const T &value)
	{
		std::ostringstream os;
		os << value;
		as_string = os.str();
	}
	const char *asString() const
	{
		return as_string.c_str();
	}
};

// points/dimensions
template <>
class ValueTraits<profit::_2dcoordinate> : public StreamedString<profit::_2dcoordinate>
{
public:
	explicit ValueTraits(const profit::_2dcoordinate &value) : StreamedString(value) {}
};

} // namespace CxxTest


namespace profit {

class LibraryInitializationFixture : public CxxTest::GlobalFixture {

private:
	std::string tmp_profit_home;

public:
	bool setUpWorld()
	{
		tmp_profit_home = ".profit_";
		while (profit::dir_exists(tmp_profit_home)) {
			tmp_profit_home += "x";
		}
		profit::setenv("PROFIT_HOME", tmp_profit_home);
		return profit::init();
	}

	bool tearDownWorld()
	{
		profit::finish();
		profit::recursive_remove(tmp_profit_home);
		return true;
	}

};

#ifndef PROFIT_TEST_NO_LIBRARY_FIXTURE
static LibraryInitializationFixture _lib_fixture;
#endif

enum zero_treatment_t {
	EXPECT_0 = 0,
	ASSUME_0
};

double relative_diff(double expected, double obtained, zero_treatment_t zero_treatment=EXPECT_0)
{
	auto diff = std::abs(expected - obtained);
	if (!diff) {
		// all good
		return 0;
	}
	// avoid NaNs due to divide-by-zero
	if (expected == 0) {
		switch (zero_treatment) {
		case EXPECT_0:
			return diff;
		case ASSUME_0:
			return 0;
		default:
			throw invalid_parameter("Invalid zero_treatment");
		}
	}
	return diff / expected;
}

void assert_images_relative_delta(const Image &expected, const Image &obtained,
    double tolerance=0, zero_treatment_t zero_treatment=EXPECT_0)
{
	// Compare dimensions, total flux, and pixel-by-pixel
	TS_ASSERT_EQUALS(expected.getDimensions(), obtained.getDimensions());
	TS_ASSERT_LESS_THAN_EQUALS(relative_diff(expected.total(), obtained.total()), tolerance);
	auto width = expected.getWidth();
	for(unsigned int i=0; i!=expected.size(); i++) {
		auto rel_diff = relative_diff(expected[i], obtained[i], zero_treatment);
		if (rel_diff > tolerance) {
			std::ostringstream msg;
			msg << "Pixel [" << i % width << "," << i / width << "] has values that are too different: ";
			msg << expected[i] << " v/s " << obtained[i];
			TS_FAIL(msg.str());
		}
	}
}

void assert_masks(const Mask &expected, const Mask &obtained)
{
	// Compare dimensions, total flux, and pixel-by-pixel
	TS_ASSERT_EQUALS(expected.getDimensions(), obtained.getDimensions());
	auto width = expected.getWidth();
	bool failed = false;
	for(unsigned int i=0; i!=expected.size(); i++) {
		std::ostringstream msg;
		msg << "Cell [" << i % width << "," << i / width << "] has different values: ";
		msg << expected[i] << " v/s " << obtained[i];
		// need explicit cast to bool here as Apple's clang refused to ultimately
		// perform a C-style cast between std::__bit_iterator and unsigned char *,
		// which is what happens under the hood in cxxtest
		TSM_ASSERT_EQUALS(msg.str(), bool(expected[i]), bool(obtained[i]));
		if (expected[i] != obtained[i]) {
			failed = true;
		}
	}
	if (failed) {
		std::ostringstream msg;
		msg << "Failed when comparing masks: Expected\n" << expected;
		msg << "\nObtained:\n" << obtained;
		TS_FAIL(msg.str());
	}
}

} // namespace profit