//
// FITS utils for profit-cli
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
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

#ifndef PROFIT_FITS_UTILS_H_
#define PROFIT_FITS_UTILS_H_

#include <string>

#include "profit/profit.h"

namespace profit {

class invalid_file : public std::exception {
private:
	std::string _what;
public:
	explicit invalid_file(const std::string &what) : _what(what) {}
	~invalid_file() noexcept {};
	const char *what() const noexcept {
		return _what.c_str();
	}
};

/// Read an image from a fits file, together with its pixel scale
Image from_fits(const std::string &filename, PixelScale &pixel_scale);

/// Export an image into a fits file, together with its pixel scale and the origin's offset
void to_fits(const Image &image, const Point &offset, const PixelScale &pixel_scale, std::string fname);

}  // namespace profit



#endif /* PROFIT_FITS_UTILS_H_ */
