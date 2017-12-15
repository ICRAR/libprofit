/**
 * Image class implementation
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

#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>

#include "profit/image.h"

namespace profit {

Mask::Mask(unsigned int width, unsigned int height) :
	surface({width, height})
{
}

Mask::Mask(Dimensions dimensions) :
	surface(dimensions)
{
}

Mask::Mask(const std::vector<bool>& data, unsigned int width, unsigned int height) :
	surface(data, {width, height})
{
}

Mask::Mask(const std::vector<bool>& data, Dimensions dimensions) :
	surface(data, dimensions)
{
}

Mask::Mask(std::vector<bool>&& data, unsigned int width, unsigned int height) :
	surface(std::move(data), {width, height})
{
}

Mask::Mask(std::vector<bool>&& data, Dimensions dimensions) :
	surface(std::move(data), dimensions)
{
}

Mask::Mask(const Mask& other) :
	surface(other)
{
}

Mask::Mask(Mask&& other) :
	surface(std::move(other))
{
}

Image::Image(unsigned int width, unsigned int height) :
	surface({width, height})
{
}

Image::Image(Dimensions dimensions) :
	surface(dimensions)
{
}

Image::Image(const std::vector<double>& data, unsigned int width, unsigned int height) :
	surface(data, {width, height})
{
}

Image::Image(const std::vector<double>& data, Dimensions dimensions) :
	surface(data, dimensions)
{
}

Image::Image(std::vector<double> &&data, unsigned int width, unsigned int height) :
	surface(std::move(data), {width, height})
{
}

Image::Image(std::vector<double>&& data, Dimensions dimensions) :
	surface(std::move(data), dimensions)
{
}

Image::Image(const Image& other) :
	surface(other)
{
}

Image::Image(Image&& other) :
	surface(std::move(other))
{
}

double Image::getTotal() const {
	const auto &data = getData();
	return std::accumulate(data.begin(), data.end(), 0.);
}

void Image::normalize()
{
	double sum = getTotal();
	if( sum > 0 ) {
		*this /= sum;
	}
}

Image Image::normalize() const
{
	Image normalized(*this);
	normalized.normalize();
	return normalized;
}

Image &Image::operator&=(const Mask &mask)
{
	// Don't apply empty masks
	if (mask.empty()) {
		return *this;
	}

	auto &data = getData();
	const auto &mask_data = mask.getData();
	std::transform(data.begin(), data.end(), mask_data.begin(), data.begin(),
		[](const double i, const bool m) {
			return m ? i : 0.;
	});
	return *this;
}

const Image Image::operator&(const Mask &mask) const
{
	Image masked(*this);
	masked &= mask;
	return masked;
}

Image &Image::operator+=(const Image& rhs)
{
	auto &data = getData();
	const auto &other_data = rhs.getData();
	std::transform(data.begin(), data.end(), other_data.begin(), data.begin(), std::plus<double>());
	return *this;
}

Image Image::operator+(const Image& rhs) const
{
	Image sum(*this);
	sum += rhs;
	return sum;

}

Image &Image::operator/=(double denominator)
{
	using std::placeholders::_1;
	auto &data = getData();
	std::transform(data.begin(), data.end(), data.begin(),
	               std::bind(std::divides<double>(), _1, denominator));
	return *this;
}

Image Image::operator/(double denominator) const
{
	Image sum(*this);
	sum /= denominator;
	return sum;
}

Image Image::operator/(int denominator) const
{
	return operator/(static_cast<double>(denominator));
}

Image Image::operator/(unsigned int denominator) const
{
	return operator/(static_cast<double>(denominator));
}

}