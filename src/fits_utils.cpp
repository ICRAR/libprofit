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

void to_fits(const Image &image, const Point &offset, const PixelScale &pixel_scale, std::string fname)
{
	FILE *f;
	char hdr[80];

	/* Append .fits if not in the name yet */
	size_t fname_size = fname.size();
	if( fname_size <= 5 || fname.compare(fname_size - 5, fname_size, ".fits") != 0 ) {
		fname = fname + ".fits";
	}

	f = fopen(fname.c_str(), "w+");
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
	fprintf(f, "%-80s", "SIMPLE  =                    T / File conforms to FITS standard");
	fprintf(f, "%-80s", "BITPIX  =                  -64 / Bits per pixel");
	fprintf(f, "%-80s", "NAXIS   =                    2 / Number of axes");
	sprintf(hdr, "NAXIS1  =           %10.0u / Width", image.getWidth());
	fprintf(f, "%-80s", hdr);
	sprintf(hdr, "NAXIS2  =           %10.0u / Height", image.getHeight());
	fprintf(f, "%-80s", hdr);
	fprintf(f, "%-80s", "CRPIX1  = 1");
	sprintf(hdr, "CRVAL1  = %f", (0.5 - offset.x) * scale_x);
	fprintf(f, "%-80s", hdr);
	sprintf(hdr, "CDELT1  = %f", scale_x);
	fprintf(f, "%-80s", hdr);
	fprintf(f, "%-80s", "CTYPE1  = ' '");
	fprintf(f, "%-80s", "CUNIT1  = ' '");
	fprintf(f, "%-80s", "CRPIX2  = 1");
	sprintf(hdr, "CRVAL2  = %f", (0.5 - offset.y) * scale_y);
	fprintf(f, "%-80s", hdr);
	sprintf(hdr, "CDELT2  = %f", scale_y);
	fprintf(f, "%-80s", hdr);
	fprintf(f, "%-80s", "CTYPE2  = ' '");
	fprintf(f, "%-80s", "CUNIT2  = ' '");
	fprintf(f, "%-80s", "END");

	auto pos = (unsigned int)ftell(f);
	auto padding = FITS_BLOCK_SIZE - (pos % FITS_BLOCK_SIZE);
	std::string spaces(padding, ' ');
	fwrite(spaces.c_str(), 1, padding, f);

	/* data has to be big-endian */
	size_t image_size = image.size();
	if( is_little_endian() ) {
		std::vector<double> big_endian_image(image_size);
		std::transform(image.begin(), image.end(), big_endian_image.begin(), swap_bytes);
		fwrite(big_endian_image.data(), sizeof(double), image_size, f);
	}
	else {
		fwrite(image.data(), sizeof(double), image_size, f);
	}

	/* Pad with zeroes until we complete the current 36*80 block */
	padding = FITS_BLOCK_SIZE - (((unsigned int)sizeof(double) * image_size) % FITS_BLOCK_SIZE);
	std::string zeros(padding, 0);
	fwrite(zeros.c_str(), 1, padding, f);
	fclose(f);
}

} // namespace profit