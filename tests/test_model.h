/**
 * model tests
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

#include "common_test_setup.h"

using namespace profit;

class TestModel : public CxxTest::TestSuite {

public:

	void test_valid_dimensions(void) {

		Model m;

		// Only the final combination is valid
		m.set_dimensions({0, 0});
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		m.set_dimensions({0, 1});
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		m.set_dimensions({1, 0});
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		m.set_dimensions({1, 1});
		m.evaluate(); // fine...

		// We generate the correct size
		Dimensions dims {100, 100};
		m.set_dimensions(dims);
		TS_ASSERT_EQUALS(dims.x * dims.y, m.evaluate().size());
	}

	void test_valid_scales(void) {

		Model m {1, 1};

		for(auto scale: {-2, -1, 0}){
			m.set_image_pixel_scale({scale, 1});
			TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		}
		m.set_image_pixel_scale({0.1, 1});
		m.evaluate(); // fine

		for(auto scale: {-2, -1, 0}){
			m.set_image_pixel_scale({0.1, scale});
			TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		}
		m.set_image_pixel_scale({0.1, 0.1});
		m.evaluate(); // fine

	}

	void test_valid_profiles(void) {

		Model m;
		for(auto name: {"unknown", "sersic1", "Sersic", "sersi", " sersic", "sersic "}) {
			TS_ASSERT_THROWS(m.add_profile(name), const invalid_parameter &);
		}
		TS_ASSERT_EQUALS(false, m.has_profiles());

		// these are fine
		for(auto name: {"brokenexp", "coresersic", "ferrer", "ferrers",
		                "king", "moffat", "psf", "sersic", "sky"}) {
			m.add_profile(name);
		}
		TS_ASSERT_EQUALS(true, m.has_profiles());

	}

	void test_valid_psf(void) {

		Model m {2, 2};

		// so far so good...
		m.evaluate();

		// add a profile, ask to convolve it, but don't give a psf
		auto skyp = m.add_profile("sky");
		skyp->parameter("convolve", true);
		TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);

		// Give a psf, everything is in order...
		m.set_psf({{1, 1}, 1, 2});
		m.evaluate();
	}

	void _add_sersic(Model &m, double xcen, double ycen, double re, bool convolve = false) {
		auto sersic = m.add_profile("sersic");
		sersic->parameter("xcen", xcen);
		sersic->parameter("ycen", ycen);
		sersic->parameter("re", re);
		sersic->parameter("convolve", convolve);
	}

	void test_profile_images_addition() {

		// three individual model images are summed up
		Model m1(100, 100);
		_add_sersic(m1, 50, 50, 10);
		auto image1 = m1.evaluate();

		Model m2(100, 100);
		_add_sersic(m2, 30, 10, 16);
		auto image2 = m2.evaluate();

		Model m3(100, 100);
		_add_sersic(m3, 23, 89, 1.2);
		auto image3 = m3.evaluate();

		// image1 holds the final result
		image1 = image1 + image2 + image3;

		// A single model image with all profile images
		Model m4(100, 100);
		_add_sersic(m4, 50, 50, 10);
		_add_sersic(m4, 30, 10, 16);
		_add_sersic(m4, 23, 89, 1.2);
		auto image4 = m4.evaluate();

		// They should be the same! We add them in the same order to make sure
		// that floating-point rounding yields the same result
		TS_ASSERT(image1 == image4);
	}

	void test_profile_images_addition_after_convolving() {

		auto convolver = create_convolver(ConvolverType::BRUTE);
		auto psf = {0., 1., 2., 3.};

		// two individual model images are summed up
		// the second is actually convolved with a psf
		Model m1(100, 100);
		_add_sersic(m1, 50, 50, 10);
		auto image1 = m1.evaluate();

		Model m2(100, 100);
		m2.set_convolver(convolver);
		m2.set_psf({psf, 2, 2});
		_add_sersic(m2, 30, 10, 16, true);
		auto image2 = m2.evaluate();

		// image1 holds the final result
		image1 = image1 + image2;

		// A single model image with all profile images
		Model m3(100, 100);
		m3.set_convolver(convolver);
		m3.set_psf({psf, 2, 2});
		_add_sersic(m3, 50, 50, 10);
		_add_sersic(m3, 30, 10, 16, true);
		auto image3 = m3.evaluate();

		// They should be the same! We add them in the same order to make sure
		// that floating-point rounding yields the same result
		TS_ASSERT(image1 == image3);

	}

	void test_finesampling()
	{

		// Finesampling produces bigger images
		Model m {100, 200};
		m.set_finesampling(2);
		TS_ASSERT_EQUALS(m.evaluate().getDimensions(), Dimensions(200, 400));
		m.set_finesampling(1);
		TS_ASSERT_EQUALS(m.evaluate().getDimensions(), Dimensions(100, 200));

		auto prepare_model = [](Model &m) {
			auto sersic = m.add_profile("sersic");
			sersic->parameter("xcen", 50.);
			sersic->parameter("ycen", 50.);
			sersic->parameter("re", 10.);
		};

		// Finesampling produces the same effect than generating a bigger model
		// image with smaller pixel scale
		Model m_orig {200, 200};
		m_orig.set_image_pixel_scale({0.5, 0.5});
		prepare_model(m_orig);
		auto im_orig = m_orig.evaluate();

		Model m_finesampled {100, 100};
		m_finesampled.set_finesampling(2);
		prepare_model(m_finesampled);
		auto im_finesampled = m_finesampled.evaluate();

		// Pixels should be equal
		auto orig_it = im_orig.begin();
		auto fine_it = im_finesampled.begin();
		for(; orig_it != im_orig.end(); orig_it++, fine_it++) {
			TS_ASSERT_EQUALS(*orig_it, *fine_it);
		}
	}

	void test_finesampling_dimensions()
	{
		Model m {100, 200};

		// no finesampling
		m.set_finesampling(1);
		TS_ASSERT_EQUALS(m.evaluate().getDimensions(), Dimensions(100, 200));

		// finesampling = 2, dimensions should be double
		m.set_finesampling(2);
		TS_ASSERT_EQUALS(m.evaluate().getDimensions(), Dimensions(200, 400));

		// finesampling = 2, return finesampled = false, dimensions should be original
		m.set_finesampling(2);
		m.set_return_finesampled(false);
		TS_ASSERT_EQUALS(m.evaluate().getDimensions(), Dimensions(100, 200));
	}

	void test_finesampling_flux()
	{
		Model m {100, 200};
		auto p = m.add_profile("sersic");
		p->parameter("xcen", 50.);
		p->parameter("ycen", 100.);
		p->parameter("re", 30.);
		auto flux = m.evaluate().total();

		// Total flux should be more or less maintained when finesampling
		m.set_finesampling(2);
		auto finesampled_flux = m.evaluate().total();
		TS_ASSERT_DELTA(flux, finesampled_flux, flux * 0.001);

		// Total flux should also be maintained when finesampling *and* not returning the fine image
		m.set_return_finesampled(false);
		finesampled_flux = m.evaluate().total();
		TS_ASSERT_DELTA(flux, finesampled_flux, flux * 0.001);
	}

	void _test_no_crop(const Dimensions &dims, unsigned int finesampling)
	{
		if (!has_fftw()) {
			TS_SKIP("No FFTW available");
		}

		auto fft_conv = create_convolver("fft", ConvolverCreationPreferences(dims * finesampling, {2, 2}, 1, nullptr, effort_t::ESTIMATE, true, simd_instruction_set::AUTO));
		Model m {dims.x, dims.y};
		m.set_convolver(fft_conv);
		m.set_psf({{1, 1, 1, 1}, 2, 2});
		m.set_finesampling(finesampling);
		auto p = m.add_profile("null");
		p->parameter("convolve", true);

		// Adjust our expectations
		auto expected_img_dims = dims * finesampling;
		auto _half_dims = expected_img_dims / 2;
		auto expected_offset = Point(_half_dims.x - 1, _half_dims.y - 1);

		// Default is to crop the image
		auto original_image = m.evaluate();
		TS_ASSERT_EQUALS(original_image.getDimensions(), expected_img_dims);

		// We know that the fft convolver creates an image that with dimensions
		// that are double of the original, and that will start at 49,49
		Point offset;
		m.set_crop(false);
		auto uncropped_image = m.evaluate(offset);
		TS_ASSERT_EQUALS(uncropped_image.getDimensions(), expected_img_dims * 2);
		TS_ASSERT_EQUALS(offset, expected_offset);

		// Manually cropping should yield the original
		TS_ASSERT_EQUALS(original_image, uncropped_image.crop(expected_img_dims, offset));
	}

	void test_no_crop()
	{
		_test_no_crop({100, 100}, 1);
	}

	void test_no_cropping_with_finesampling()
	{
		_test_no_crop({100, 100}, 2);
	}

};