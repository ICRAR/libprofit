/**
 * Image class tests
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2017
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

#include "profit/profit.h"

using namespace profit;

class TestImage : public CxxTest::TestSuite {

private:
	void assert_empty(const Image &im) {
		TS_ASSERT_EQUALS(true, im.empty());
		TS_ASSERT_EQUALS(0, im.getData().size());
		TS_ASSERT_EQUALS(0, im.getWidth());
		TS_ASSERT_EQUALS(0, im.getHeight());
		TS_ASSERT_EQUALS(0, im.getSize());
	}

	void assert_not_empty(const Image &im) {
		TS_ASSERT_DIFFERS(true, im.empty());
		TS_ASSERT_DIFFERS(0, im.getData().size());
		TS_ASSERT_DIFFERS(0, im.getWidth());
		TS_ASSERT_DIFFERS(0, im.getHeight());
		TS_ASSERT_DIFFERS(0, im.getSize());
	}

public:

	void test_empty() {
		Image empty;
		assert_empty(empty);
	}

	void test_correct_dimensions() {
		Image im1;
		Image im2({1}, 1, 1);
		Image im3({1, 2}, 1, 2);
		Image im4({1, 2}, 2, 1);
		Image im5({1, 2, 3}, 1, 3);
		Image im6({1, 2, 3}, 3, 1);
		Image im7({1, 2, 3, 4}, 2, 2);
	}

	void test_invalid_dimensions() {
		TS_ASSERT_THROWS(Image im({1, 2, 3, 4}, 1, 1), std::invalid_argument);
	}

	void test_copy() {

		// Copy assignment
		Image im1({1,2,3,4}, 2, 2);
		assert_not_empty(im1);
		Image im2 = im1;
		assert_not_empty(im1);
		assert_not_empty(im2);
		TS_ASSERT(im1 == im2);

		// Copy ctor
		Image im3({1,2,3,4}, 2, 2);
		assert_not_empty(im3);
		Image im4(im3);
		assert_not_empty(im3);
		assert_not_empty(im4);
		TS_ASSERT(im3 == im4);

		// Construct from vector with data copying
		std::vector<double> data {1, 2, 3, 4};
		Image im5(data, 2, 2);
		assert_not_empty(im5);
		TS_ASSERT_EQUALS(false, data.empty());
		TS_ASSERT_EQUALS(data, im5.getData());
	}

	void test_move() {

		// Move assignment
		Image im1({1,2,3,4}, 2, 2);
		assert_not_empty(im1);
		Image im2(std::move(im1));
		assert_empty(im1);
		assert_not_empty(im2);

		// Move ctor
		Image im3({1,2,3,4}, 2, 2);
		assert_not_empty(im3);
		Image im4 = std::move(im3);
		assert_empty(im3);
		assert_not_empty(im4);

		// Construct from vector with data movement
		std::vector<double> data {1, 2, 3, 4};
		Image im5(std::move(data), 2, 2);
		assert_not_empty(im5);
		TS_ASSERT_EQUALS(true, data.empty());
	}

	void test_add() {

		std::vector<double> data {1, 2, 3, 4};
		Image im1(data, 2, 2);
		Image im2(data, 2, 2);

		// Normal sum
		Image im3 = im1 + im2;
		// += sum
		im2 += im1;

		// Both sums should have worked
		for(auto im: {&im3, &im2}) {
			for(unsigned int i = 0; i < 4; i++) {
				TS_ASSERT_DELTA(im->getData()[i], data[i] * 2, 1e-6);
			}
		}
	}

	void test_divide() {

		std::vector<double> data {1, 2, 3, 4};
		Image im1(data, 2, 2);

		Image im2 = im1 / 2;
		im1 /= 2;

		// Both divisions should have worked
		for(auto im: {&im1, &im2}) {
			for(unsigned int i = 0; i < 4; i++) {
				TS_ASSERT_DELTA(im->getData()[i], data[i] / 2, 1e-6);
			}
		}
	}

	void test_normalize() {

		// use both flavours: const and not const
		Image im1({1, 2, 3, 4}, 2, 2);
		Image im2 = static_cast<const Image &>(im1).normalize();
		im1.normalize();
		for(auto im: {&im1, &im2}) {
			TS_ASSERT_DELTA(1, im->getTotal(), 1e-6);
		}

		// A zero-values image doesn't get normalized
		Image im3({0, 0, 0, 0}, 2, 2);
		Image im4 = static_cast<const Image &>(im3).normalize();
		im3.normalize();
		for(auto im: {&im3, &im4}) {
			TS_ASSERT_DELTA(0, im->getTotal(), 1e-6);
		}

		// Re-normalizing should get us to the same place
		const Image im5({1, 2, 3, 4}, 2, 2);
		Image im6 = static_cast<const Image &>(im5.normalize()).normalize();
		TS_ASSERT_DELTA(1, im6.getTotal(), 1e-6);
	}

	void test_crop() {


		// Dimensions are correct
		Image im1({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 5, 2);
		Image im2 = im1.crop(2, 2, 0, 0);
		TS_ASSERT_EQUALS(2, im2.getWidth());
		TS_ASSERT_EQUALS(2, im2.getHeight());
		im2 = im1.crop(4, 1, 0, 1);
		TS_ASSERT_EQUALS(4, im2.getWidth());
		TS_ASSERT_EQUALS(1, im2.getHeight());

		// Try all combinations now
		im1 = Image({1, 2, 3, 4}, 2, 2);

		// Should get one pixel at a time
		im2 = im1.crop(1, 1, 0, 0);
		TS_ASSERT_DELTA(1, im2.getData()[0], 1e-6);
		im2 = im1.crop(1, 1, 1, 0);
		TS_ASSERT_DELTA(2, im2.getData()[0], 1e-6);
		im2 = im1.crop(1, 1, 0, 1);
		TS_ASSERT_DELTA(3, im2.getData()[0], 1e-6);
		im2 = im1.crop(1, 1, 1, 1);
		TS_ASSERT_DELTA(4, im2.getData()[0], 1e-6);

		// should get the one row at a time
		im2 = im1.crop(2, 1, 0, 0);
		TS_ASSERT_DELTA(1, im2.getData()[0], 1e-6);
		TS_ASSERT_DELTA(2, im2.getData()[1], 1e-6);
		im2 = im1.crop(2, 1, 0, 1);
		TS_ASSERT_DELTA(3, im2.getData()[0], 1e-6);
		TS_ASSERT_DELTA(4, im2.getData()[1], 1e-6);

		// should get one column at a time
		im2 = im1.crop(1, 2, 0, 0);
		TS_ASSERT_DELTA(1, im2.getData()[0], 1e-6);
		TS_ASSERT_DELTA(3, im2.getData()[1], 1e-6);
		im2 = im1.crop(1, 2, 1, 0);
		TS_ASSERT_DELTA(2, im2.getData()[0], 1e-6);
		TS_ASSERT_DELTA(4, im2.getData()[1], 1e-6);

		// gotta get them all!
		im2 = im1.crop(2, 2, 0, 0);
		TS_ASSERT(im2 == im1);
	}

	void test_invalid_crops() {

		Image im({1, 2, 3, 4}, 2, 2);

		// too wide
		TS_ASSERT_THROWS(im.crop(3, 1, 0, 0), std::invalid_argument);

		// too tall
		TS_ASSERT_THROWS(im.crop(1, 3, 0, 0), std::invalid_argument);

		// horizontally overflows
		TS_ASSERT_THROWS(im.crop(0, 0, 3, 1), std::invalid_argument);

		// vertically overflows
		TS_ASSERT_THROWS(im.crop(0, 0, 1, 3), std::invalid_argument);
	}

	void test_extend() {

		Image im({1, 2, 3, 4}, 2, 2);

		// Correct dimensions
		for(auto x: {4, 10}) {
			for(auto y: {4, 10}) {
				Image im2 = im.extend(x, y, 0, 0);
				TS_ASSERT_EQUALS(im2.getWidth(), x);
				TS_ASSERT_EQUALS(im2.getHeight(), y);
			}
		}

		// leave at 0, 0
		auto im2 = im.extend(3, 3, 0, 0);
		TS_ASSERT_EQUALS(im.getData()[0], im2.getData()[0]);
		TS_ASSERT_EQUALS(im.getData()[1], im2.getData()[1]);
		TS_ASSERT_EQUALS(im.getData()[2], im2.getData()[3]);
		TS_ASSERT_EQUALS(im.getData()[3], im2.getData()[4]);

		// leave at 1, 1
		im2 = im.extend(3, 3, 1, 1);
		TS_ASSERT_EQUALS(im.getData()[0], im2.getData()[4]);
		TS_ASSERT_EQUALS(im.getData()[1], im2.getData()[5]);
		TS_ASSERT_EQUALS(im.getData()[2], im2.getData()[7]);
		TS_ASSERT_EQUALS(im.getData()[3], im2.getData()[8]);

		// extend to the same dimensions, should be equals
		im2 = im.extend(2, 2, 0, 0);
		TS_ASSERT(im == im2);
	}

	void test_invalid_extends() {

		Image im({1, 2, 3, 4}, 2, 2);

		// too thin
		TS_ASSERT_THROWS(im.extend(1, 5, 0, 0), std::invalid_argument);

		// too short
		TS_ASSERT_THROWS(im.crop(5, 1, 0, 0), std::invalid_argument);

		// horizontally overflows
		TS_ASSERT_THROWS(im.crop(5, 5, 4, 0), std::invalid_argument);

		// vertically overflows
		TS_ASSERT_THROWS(im.crop(5, 5, 0, 4), std::invalid_argument);

	}
};