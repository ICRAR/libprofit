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

} // namespace profit