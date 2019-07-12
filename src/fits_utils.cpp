//
// FITS utilities for profit-cli
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

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <utility>

#include "profit/fits_utils.h"
#include "profit/model.h"

namespace profit {

static inline
bool is_little_endian() {
	const uint32_t i = 0x01234567;
	return *reinterpret_cast<const uint8_t *>(&i) == 0x67;
}

static inline
double swap_bytes(const double v) {
	double r;
	const char *vbytes = reinterpret_cast<const char *>(&v);
	char *rbytes = reinterpret_cast<char *>(&r);
	rbytes[0] = vbytes[7];
	rbytes[1] = vbytes[6];
	rbytes[2] = vbytes[5];
	rbytes[3] = vbytes[4];
	rbytes[4] = vbytes[3];
	rbytes[5] = vbytes[2];
	rbytes[6] = vbytes[1];
	rbytes[7] = vbytes[0];
	return r;
}

static constexpr unsigned int FITS_BLOCK_SIZE = 36 * 80;

Image from_fits(const std::string &filename, PixelScale &pixel_scale)
{
	std::ifstream f(filename, std::ios_base::binary);
	if (!f) {
		std::ostringstream os;
		os << "Couldn't open '" << filename << "' for reading: " << std::strerror(errno);
		throw invalid_file(os.str());
	}

	/*
	 * Standard headers, we're assuming they say 'T" for SIMPLE, -64 for BITPIX
	 * and 2 for NAXIS.
	 */
	char hdr[80];
	unsigned int width = 0;
	unsigned int height = 0;
	double scale_x = 1;
	double scale_y = 1;
	while (true) {

		f.read(hdr, 80);
		if (f.gcount() != 80) {
			std::ostringstream os;
			os << "File " << filename << " does not look like a normal fits file";
			throw invalid_file(os.str());
		}

		if( !std::strncmp("NAXIS1", hdr, 6) ) {
			std::sscanf(hdr, "NAXIS1 = %u", &width);
		}
		else if( !std::strncmp("NAXIS2", hdr, 6) ) {
			std::sscanf(hdr, "NAXIS2 = %u", &height);
		}
		else if( !std::strncmp("CDELT1", hdr, 6) ) {
			std::sscanf(hdr, "CDELT1 = %lf", &scale_x);
		}
		else if( !std::strncmp("CDELT2", hdr, 6) ) {
			std::sscanf(hdr, "CDELT2 = %lf", &scale_y);
		}
		else if( !std::strncmp("END", hdr, 3) ) {
			break;
		}
	}

	pixel_scale = {scale_x, scale_y};

	// Move until the end of the header and read the actual data
	auto pos = f.tellg();
	auto padding = FITS_BLOCK_SIZE - (pos % FITS_BLOCK_SIZE);
	f.seekg(padding, std::ios_base::cur);

	Image psf(width, height);
	std::size_t n_bytes = psf.size() * sizeof(double);
	f.read(reinterpret_cast<char *>(psf.data()), n_bytes);
	if (std::size_t(f.gcount()) != n_bytes) {
		throw invalid_file("Error while reading file: less data found than expected");
	}
	f.close();

	/* data has to be big-endian */
	if( is_little_endian() ) {
		std::transform(psf.begin(), psf.end(), psf.begin(), swap_bytes);
	}

	return psf;
}

static
void write_header(std::ofstream &f, const char *header) {
	char hdr[81];
	std::sprintf(hdr, "%-80s", header);
	f.write(hdr, 80);
};

template <typename ...Ts>
static
void write_header(std::ofstream &f, const char *fmt, Ts&&...values) {
	char hdr[81];
	std::sprintf(hdr, fmt, std::forward<Ts>(values)...);
	write_header(f, hdr);
};


void add_padding(std::ofstream &f, std::ofstream::pos_type from, char padding_char)
{
	auto pad_size = FITS_BLOCK_SIZE - (from % FITS_BLOCK_SIZE);
	std::string padding(pad_size, padding_char);
	f.write(padding.c_str(), padding.size());
}

void to_fits(const Image &image, const Point &offset, const PixelScale &pixel_scale, std::string fname)
{
	/* Append .fits if not in the name yet */
	size_t fname_size = fname.size();
	if( fname_size <= 5 || fname.compare(fname_size - 5, fname_size, ".fits") != 0 ) {
		fname = fname + ".fits";
	}

	std::ofstream f(fname, std::ios_base::binary);
	if (!f) {
		std::ostringstream os;
		os << "Couldn't open '" << fname << "' for writing: " << std::strerror(errno);
		throw invalid_file(os.str());
	}

	/*
	 * Standard headers
	 *
	 * The first five headers are required, and must be in "fixed format",
	 * meaning that their values must be right-indented on column 30, sigh...
	 */
	auto scale_x = pixel_scale.first;
	auto scale_y = pixel_scale.second;
	write_header(f, "SIMPLE  =                    T / File conforms to FITS standard");
	write_header(f, "BITPIX  =                  -64 / Bits per pixel");
	write_header(f, "NAXIS   =                    2 / Number of axes");
	write_header(f, "NAXIS1  =           %10.0u / Width", image.getWidth());
	write_header(f, "NAXIS2  =           %10.0u / Height", image.getHeight());
	write_header(f, "CRPIX1  = 1");
	write_header(f, "CRVAL1  = %f", (0.5 - offset.x) * scale_x);
	write_header(f, "CDELT1  = %f", scale_x);
	write_header(f, "CTYPE1  = ' '");
	write_header(f, "CUNIT1  = ' '");
	write_header(f, "CRPIX2  = 1");
	write_header(f, "CRVAL2  = %f", (0.5 - offset.y) * scale_y);
	write_header(f, "CDELT2  = %f", scale_y);
	write_header(f, "CTYPE2  = ' '");
	write_header(f, "CUNIT2  = ' '");
	write_header(f, "END");
	add_padding(f, f.tellp(), ' ');

	/* data has to be big-endian */
	size_t image_size = image.size();
	if (is_little_endian()) {
		std::for_each(image.begin(), image.end(), [&f](double pixel) {
			pixel = swap_bytes(pixel);
			f.write(reinterpret_cast<const char *>(&pixel), sizeof(double));
		});
	}
	else {
		f.write(reinterpret_cast<const char *>(image.data()), sizeof(double) * image_size);
	}

	add_padding(f, sizeof(double) * image_size, 0);
}

} // namespace profit