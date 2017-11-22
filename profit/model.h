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
#include <utility>
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
 * A pair with an Image in the first element, and a Point representing an
 * offset in the Image in the second element.
 */
typedef std::pair<Image, Point> ImageAndOffset;

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
	 * The result of the computation is returned as an Image (the first element
	 * of the pair returned by this function), which may be of a different size
	 * from the one originally requested if the user set this model's @ref crop
	 * property to ``false``. The second element corresponds to the offset
	 * at which the image resulting of evaluating this Model with its configured
	 * parameters is with respect to the Image value returned in the first
	 * element.
	 *
	 * In other words, the Image in the first element can be bigger than the
	 * Model's @ref width and @ref height if the user requested this Model to
	 * return a non-cropped image.
	 *
	 * @returns A pair with the image created by libprofit in as first element
	 *          and an offset in the second.
	 */
	ImageAndOffset evaluate();

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
	ConvolverPtr convolver;

	/**
	 * Due to their internal workings, some convolvers produce actually bigger
	 * which are (by default) cropped to the size of the original images created
	 * by the profiles. If this option is set to true, then the result of the
	 * convolution will *not* be cropped, meaning that the result of the model
	 * evaluation will be bigger than what was originally requested.
	 */
	bool crop;

	/**
	 * Whether the actual evaluation of profiles should be skipped or not.
	 * Profile validation still occurs.
	 */
	bool dry_run;

#ifdef PROFIT_OPENCL
	OpenCLEnvPtr opencl_env;
#endif /* PROFIT_OPENCL */

#ifdef PROFIT_OPENMP
	/**
	 * Maximum number of OpenMP threads to use to evaluate the profiles
	 * contained in this model. 0 threads means that no OpenMP support
	 * has been requested.
	 */
	unsigned int omp_threads;
#endif /* PROFIT_OPENMP */

private:

	/**
	 * A list of pointers to the individual profiles used to generate the
	 * model's image
	 */
	std::vector<std::shared_ptr<Profile>> profiles;

};

} /* namespace profit */

#endif /* PROFIT_MODEL_H */
