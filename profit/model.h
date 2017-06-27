/**
 * Header file for the main Model class of libprofit
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

#ifndef PROFIT_MODEL_H
#define PROFIT_MODEL_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "profit/config.h"
#include "profit/convolve.h"
#include "profit/fft.h"
#include "profit/opencl.h"

namespace profit
{

/* Forward declaration */
struct ProfileStats;
class Profile;

/**
 * The overall model to be created
 *
 * The model includes the width and height of the image to produce, as well as
 * the resolution to use when performing calculations. Having resolution
 * allows us to specify pixel position with decimal places; e.g., the center
 * point for a given profile.
 */
class Model {

public:

	/**
	 * The types of convolvers this Model can internally create
	 */
	enum ConvolverType {
		BRUTE = 0,
#ifdef PROFIT_OPENCL
		OPENCL,
#endif // PROFIT_OPENCL
#ifdef PROFIT_FFTW
		FFT,
#endif // PROFIT_FFTW
	};

	/**
	 * Constructor
	 *
	 * It creates a new model to which profiles can be added, and that can be
	 * used to calculate an image.
	 */
	Model(unsigned int width = 0, unsigned int height = 0);

	/**
	 * Creates a new profile for the given name and adds it to the given model.
	 * On success, the new profile is created, added to the model,
	 * and its reference is returned for further customization.
	 * If a profile with the given name is not supported an invalid_parameter
	 * exception is thrown.
	 *
	 * @param profile_name The name of the profile that should be created
	 * @returns A shared pointer to the new profile that corresponds to the given name
	 */
	std::shared_ptr<Profile> add_profile(const std::string &profile_name);

	/**
	 * Whether this model contains any profiles or not.
	 *
	 * @return `true` if this module contains at least one profile,
	 * `false` otherwise
	 */
	bool has_profiles() const;

	/**
	 * Calculates an image using the information contained in the model.
	 * The result of the computation is stored in the image field.
	 *
	 * @returns The image created by libprofit. The data is organized by rows
	 *          first, columns later; i.e pixel ``(x,y)`` is accessed by
	 *          ``image[y*width + x]``
	 */
	std::vector<double> evaluate();

	/**
	 * Creates and returns a new Convolver instance suitable to be used by this
	 * Model's settings (i.e., image width/height, psf width/heigh, and number
	 * of OpenMP threads, if supported).
	 *
	 * This method does not set the internal convolver to the returned instance;
	 * this must be done by the user.
	 *
	 * @return A convolver suitable to be used with this model's settings.
	 */
	std::shared_ptr<Convolver> create_convolver() const;

#ifdef PROFIT_DEBUG
	std::map<std::string, std::map<int, int>> get_profile_integrations() const;
#endif

	/**
	 * Return a map of all profile statistics.
	 *
	 * @return A map indexed by profile name with runtime statistics
	 */
	std::map<std::string, std::shared_ptr<ProfileStats>> get_stats() const;

	/**
	 * The width of the model to generate
	 */
	unsigned int width;

	/**
	 * The height of the model to generate
	 */
	unsigned int height;

	/**
	 * The X scale; that is, the width of a single pixel in image coordinates
	 */
	double scale_x;

	/**
	 * The Y scale; that is, the height of a single pixel in image coordinates
	 */
	double scale_y;

	/**
	 * The base magnitude applied to all models
	 */
	double magzero;

	/**
	 * The point spread function (psf) to use when convolving images
	 */
	std::vector<double> psf;

	/**
	 * The psf's width
	 */
	unsigned int psf_width;

	/**
	 * The psf's height
	 */
	unsigned int psf_height;

	/**
	 * The PSF's X scale; that is, the width of a single PSF pixel in image
	 * coordinates
	 */
	double psf_scale_x;

	/**
	 * The PSF's Y scale; that is, the height of a single PSF pixel in image
	 * coordinates
	 */
	double psf_scale_y;

	/**
	 * The calculation mask. If given it must be the same size of the expected
	 * output image, and its values are used to limit the profile calculation
	 * only to a given area (i.e., those cells where the value is ``true``).
	 */
	std::vector<bool> calcmask;

	/**
	 * The object used to carry out the convolution, if necessary.
	 * If a convolver is present before calling `evaluate` then it is used.
	 * If missing, then a new one is created internally.
	 */
	std::shared_ptr<Convolver> convolver;

	/**
	 * Which type of convolver should be constructed if one is required,
	 * but missing.
	 */
	ConvolverType convolver_type;

	/**
	 * Whether the actual evaluation of profiles should be skipped or not.
	 * Profile validation still occurs.
	 */
	bool dry_run;

#ifdef PROFIT_OPENCL
	std::shared_ptr<OpenCL_env> opencl_env;
#endif /* PROFIT_OPENCL */

#ifdef PROFIT_OPENMP
	/**
	 * Maximum number of OpenMP threads to use to evaluate the profiles
	 * contained in this model. 0 threads means that no OpenMP support
	 * has been requested.
	 */
	unsigned int omp_threads;
#endif /* PROFIT_OPENMP */

#ifdef PROFIT_FFTW
	/**
	 * Whether or not a copy of the FFT'd version of the PSF should be kept
	 * by the FFT-based convolver (if used).
	 */
	bool reuse_psf_fft;

	/**
	 * How much effort should be used to create the FFT plans used by the
	 * convolver. Used only if `use_fft` is set and no convolver is set.
	 */
	FFTPlan::effort_t fft_effort;
#endif /* PROFIT_FFTW */

private:

	/**
	 * A list of pointers to the individual profiles used to generate the
	 * model's image
	 */
	std::vector<std::shared_ptr<Profile>> profiles;

};

} /* namespace profit */

#endif /* PROFIT_MODEL_H */
