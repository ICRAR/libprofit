#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2019
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of libprofit.
#
# libprofit is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# libprofit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with libprofit.  If not, see <http://www.gnu.org/licenses/>.

add_executable(profit-cli fits_utils.cpp profit-cli.cpp)
target_link_libraries(profit-cli profit GSL::gsl)
if (MSVC)
	target_include_directories(profit-cli PRIVATE ${PROFIT_GETOPT_DIR})
endif()

install(TARGETS profit-cli
        EXPORT profit-targets
        RUNTIME DESTINATION bin
)

if (LIBPROFIT_TEST)
	add_test("cli-version" profit-cli -V)
	add_test("cli-help" profit-cli -h)
	add_test("cli-null-profile" profit-cli -p null)
	add_test("cli-many-profiles" profit-cli -p 100,null -p 5,sky)
	add_test("cli-unknown-parameter" profit-cli -p null:unknown=1)
	add_test("cli-text-output" profit-cli -p null -t)
	add_test("cli-fits-output" profit-cli -p null -f image)
	add_test("cli-convolution-with-cmdline" profit-cli -p null:convolve=1 -P 2:2:1,1,2,2)
	add_test("cli-complex-usage" profit-cli -p null -w 100 -H 50 -S 2 -F -m 1 -i 2 -u -T brute -e 1 -r)
	# Needs 'cli-fits-output' to be run first
	add_test("cli-convolution-with-file" profit-cli -p null:convolve=1 -P image.fits)
	if (PROFIT_OPENCL)
		add_test("cli-opencl-list" profit-cli -c)
	endif()
endif()
