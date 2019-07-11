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

	void test_valid_masks(void)
	{
		Dimensions image_dims{100, 100};
		Dimensions psf_dims{4, 4};
		Model m(image_dims);
		auto add_convolution = [&m, &psf_dims]() {
			auto profile = m.add_profile("null");
			profile->parameter("convolve", true);
			m.set_psf(Image{1, psf_dims});
		};
		auto common_checks = [&m, &image_dims]() {
			m.set_mask(Mask{true, image_dims});
			m.evaluate(); // fine...
			m.set_mask(Mask{});
			m.evaluate(); // fine...
			m.set_mask(Mask{true, image_dims - 1});
			TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		};
		auto convolution_checks = [&m, &image_dims, &psf_dims]() {
			// pre-computed masks should have psf_padding if necessary.
			// Pre-adjusted mask loose the original coverage information,
			// so the model assumes its dimensionality reflects whether the
			// model should grow or not
			m.set_adjust_mask(false);
			m.set_mask(Mask{1, image_dims});
			m.evaluate(); // fine...
			m.set_mask(Mask{1, image_dims + psf_dims});
			m.evaluate(); // fine...
			m.set_mask(Mask{1, image_dims + psf_dims + 1});
			TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
			m.set_mask(Mask{1, image_dims + psf_dims - 1});
			TS_ASSERT_THROWS(m.evaluate(), const invalid_parameter &);
		};

		common_checks();
		add_convolution();
		common_checks();
		convolution_checks();
	}

	void test_precomputed_mask_works()
	{
		// Two images are generated: one with the pre-computed mask
		// and another with the interally-computed mask.
		// They should be the same
		Dimensions image_dims{10, 10};
		Dimensions psf_dims{4, 4};
		Image psf{1, psf_dims};
		Model m{image_dims};
		m.set_psf(psf);
		auto sersic = m.add_profile("sersic");
		sersic->parameter("xcen", 5.);
		sersic->parameter("ycen", 5.);
		sersic->parameter("re", 3.);
		sersic->parameter("convolve", true);

		auto assert_precomputed_mask_works = [&m, &image_dims, &psf](const Mask &mask) {
			m.set_adjust_mask(true);
			m.set_mask(mask);
			auto internal_mask_image = m.evaluate();

			Mask adjusted_mask(mask);
			Model::adjust(adjusted_mask, image_dims, psf);
			m.set_adjust_mask(false);
			m.set_mask(adjusted_mask);
			auto offline_mask_image = m.evaluate();

			assert_images_relative_delta(internal_mask_image, offline_mask_image & mask, 0, zero_treatment_t::EXPECT_0);
		};

		// A full mask that requires the model and mask to expand their dimensions
		assert_precomputed_mask_works(Mask{true, image_dims});

		// A small mask that doesn't require model and mask dimensions expansion,
		// even when the mask's coverage is expanded by psf_dims/2
		Mask mask1{true, image_dims - psf_dims * 2};
		mask1 = mask1.extend(image_dims, (image_dims - mask1.getDimensions()) / 2);
		assert_precomputed_mask_works(mask1);

		// A mask that doesn't require model and mask dimensions expansion
		// due to its original coverage, but would require if the Model
		// tried to add padding *again* due to the already-expanded coverage
		Mask mask2{true, image_dims - psf_dims};
		mask2 = mask2.extend(image_dims, (image_dims - mask2.getDimensions()) / 2);
		assert_precomputed_mask_works(mask2);
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
		assert_images_relative_delta(image1, image4, 1e-9, zero_treatment_t::ASSUME_0);
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

	void _test_no_crop(const Dimensions &dims, unsigned int finesampling, ConvolverType conv_type)
	{
		ConvolverCreationPreferences conv_prefs;
		conv_prefs.src_dims = dims * finesampling;
		conv_prefs.krn_dims = {2, 2};
		conv_prefs.reuse_krn_fft = true;
		Model m {dims.x, dims.y};
		m.set_convolver(create_convolver(conv_type, conv_prefs));
		m.set_psf({{1, 1, 1, 1}, 2, 2});
		m.set_finesampling(finesampling);
		auto p = m.add_profile("null");
		p->parameter("convolve", true);


		// Default is to crop the image
		auto expected_img_dims = dims * finesampling;
		auto original_image = m.evaluate();
		TS_ASSERT_EQUALS(original_image.getDimensions(), expected_img_dims);

		// Don't crop now, we should be able to reproduce that crop
		auto expected_offset = Dimensions{0, 0};
		auto expected_uncropped_img_dims = expected_img_dims;
		if (conv_type == ConvolverType::FFT) {
			expected_offset = expected_img_dims / 2 - 1;
			expected_uncropped_img_dims = expected_img_dims * 2;
		}
		Point offset;
		m.set_crop(false);
		auto uncropped_image = m.evaluate(offset);
		TS_ASSERT_EQUALS(uncropped_image.getDimensions(), expected_uncropped_img_dims);
		TS_ASSERT_EQUALS(offset, expected_offset);
		TS_ASSERT_EQUALS(original_image, uncropped_image.crop(expected_img_dims, offset));
	}

	void test_no_crop()
	{
		_test_no_crop({20, 20}, 1, ConvolverType::BRUTE);
	}

	void test_no_cropping_with_finesampling()
	{
		_test_no_crop({20, 20}, 2, ConvolverType::BRUTE);
	}

	void test_no_crop_fft_convolver()
	{
		if (!has_fftw()) {
			TS_SKIP("No FFTW available");
		}
		_test_no_crop({20, 20}, 1, ConvolverType::FFT);
	}

	void test_no_cropping_with_finesampling_fft_convolver()
	{
		if (!has_fftw()) {
			TS_SKIP("No FFTW available");
		}
		_test_no_crop({20, 20}, 2, ConvolverType::FFT);
	}

};

class TestMaskAdjustments : public CxxTest::TestSuite {

public:

	void test_null_mask()
	{
		Mask mask;
		Model::adjust(mask, {100, 100}, Image{});
		TS_ASSERT(!mask);
	}

	void test_no_dimension_extension()
	{
		Dimensions image_dims{5, 5};
		Dimensions psf_dims{2, 2};

		// A mask with the outer part unset, but small enough for it to not
		// require dimension extension
		Mask mask{true, image_dims - psf_dims};
		mask = mask.extend(image_dims, psf_dims / 2);

		Image psf{1, psf_dims};
		Model::adjust(mask, image_dims, psf);
		assert_masks(Mask{true, image_dims}, mask);
	}

	void test_with_dimension_extension()
	{
		Dimensions image_dims{5, 5};
		Dimensions psf_dims{2, 2};

		// A mask fully set, it will require dimension extensions
		Mask mask{true, image_dims};
		Image psf{1, psf_dims};
		Model::adjust(mask, image_dims, psf);
		assert_masks(Mask{true, image_dims + psf_dims}, mask);
	}

};

class TestFluxCapturing : public CxxTest::TestSuite {

private:

	void _test_flux_is_captured(unsigned int model_size, unsigned int mask_size,
	    unsigned int psf_size, ConvolverPtr convolver, unsigned int finesampling)
	{
		// Make a model containing a sky profile with a fixed background. Then
		// we evaluate the model and then we do the same but with a mask.
		// Optionally convolution with a psf is involved in both steps.
		// At the end we should end up with the same image within the masked
		// region, which implies we correctly captured all flux, even when doing
		// convolution.
		//
		// We repeat the above for all the possible different positions of the
		// mask to make sure we are accounting for all possible edge cases.
		// TODO: add for loops for crop/no-crop

		const Mask original_mask{true, mask_size, mask_size};
		Image psf{1, psf_size, psf_size};
		Model m{model_size, model_size};
		m.set_finesampling(finesampling);
		m.set_return_finesampled(false);
		auto sersic = m.add_profile("sky");
		sersic->parameter("bg", 1.);
		if (convolver) {
			sersic->parameter("convolve", true);
		}
		m.set_psf(psf);
		m.set_convolver(convolver);
		auto non_masked_image = m.evaluate();

		auto assert_same_masked_image = [&] (Point mask_offset) {
			auto mask = original_mask.extend({model_size, model_size}, {mask_offset.x, mask_offset.y});
			m.set_mask(mask);
			auto masked_image = m.evaluate();
			assert_images_relative_delta(non_masked_image & mask, masked_image, 0.0001, zero_treatment_t::EXPECT_0);
		};

		for (auto mask_y_offset = 0U; mask_y_offset <= model_size - mask_size; mask_y_offset++) {
			for (auto mask_x_offset = 0U; mask_x_offset <= model_size - mask_size; mask_x_offset++) {
				assert_same_masked_image({mask_x_offset, mask_y_offset});
			}
		}
	}

	void _test_flux_is_captured(unsigned int model_size, unsigned int mask_size, unsigned int psf_size)
	{
		std::vector<ConvolverPtr> convolvers;
		convolvers.emplace_back(nullptr);
		convolvers.emplace_back(create_convolver(ConvolverType::BRUTE));
		if (has_fftw()) {
			convolvers.emplace_back(create_convolver(ConvolverType::FFT));
		}
		if (has_opencl() && !get_opencl_info().empty()) {
			ConvolverCreationPreferences prefs;
			prefs.opencl_env = get_opencl_environment(0, 0, false, false);
			convolvers.emplace_back(create_convolver(ConvolverType::OPENCL, prefs));
		}
		for (auto &convolver: convolvers) {
			for (auto finesampling: {1, 2}) {
				_test_flux_is_captured(model_size, mask_size, psf_size, convolver, finesampling);
			}
		}
	}

public:
	void test_empty_mask()
	{
		_test_flux_is_captured(20, 0, 5);
	}

	void test_tiny()
	{
		_test_flux_is_captured(5, 2, 5);
		_test_flux_is_captured(5, 3, 5);
		_test_flux_is_captured(5, 4, 5);
		_test_flux_is_captured(5, 3, 4);
		_test_flux_is_captured(5, 4, 4);
		_test_flux_is_captured(5, 5, 4);
		_test_flux_is_captured(3, 1, 2);
	}

	void test_small_mask_small_psf()
	{
		_test_flux_is_captured(20, 4, 5);
		_test_flux_is_captured(20, 5, 5);
		_test_flux_is_captured(20, 6, 5);
		_test_flux_is_captured(20, 5, 6);
		_test_flux_is_captured(20, 6, 6);
		_test_flux_is_captured(20, 7, 6);
	}

	void test_small_mask_medium_psf()
	{
		_test_flux_is_captured(20, 5, 9);
	}

	void test_small_mask_full_psf()
	{
		_test_flux_is_captured(20, 5, 20);
	}

	void test_small_mask_bigger_psf()
	{
		_test_flux_is_captured(20, 5, 30);
	}

};