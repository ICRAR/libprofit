//
// FFT-related class definitions
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2017
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

#ifndef PROFIT_FFT_IMPL_H
#define PROFIT_FFT_IMPL_H

#include "profit/config.h"

#ifdef PROFIT_FFTW

#include <complex>
#include <memory>
#include <vector>

#include <fftw3.h>

#include "profit/common.h"
#include "profit/fft.h"

namespace profit {

/**
 * A deleter class that we can use for our unique_ptr objects holding
 * FFTW-allocated arrays.
 */
template <typename T>
class fftw_deleter {
public:
	void operator()(T *x) {
		fftw_free(x);
	}
};

/**
 * An FFT transformer that turns real numbers into complex numbers and back.
 *
 * Instances of this class are able to perform forward FFT transformation
 * from a collection of double values into its corresponding complex series and
 * back using hermitian redundancy. The forward and backward transformations
 * occur in an internal buffer, but they allow users to pass collections of the
 * corresponding type to store the output of the operation, thus avoiding
 * dynamic memory allocation.
 */
class FFTRealTransformer {

public:

	/**
	 * Creates a new transformer that will work with images and vectors of size
	 * @p size, using effort @p effort and @p omp_threads threads.
	 *
	 * @param size The size of the data to be transformed
	 * @param effort The kind of effort that should be put into creating this plan
	 * @param omp_threads The number of threads to use to execute the plan
	 */
	FFTRealTransformer(unsigned int size, effort_t effort, unsigned int omp_threads);

	/**
	 * Destructor. It destroys the underlying plans.
	 */
	~FFTRealTransformer();

	/**
	 * Transforms a container of numbers into their Fourier Transform. The
	 * resulting vector is a vector of complex values.
	 *
	 * @param input The numbers to transform
	 * @param output The vector of complex numbers where the result will be stored.
	 * It must be at least as long as the hermitian size
	 */
	template <typename T>
	void forward(const T &input, std::vector<std::complex<double>> &output) const;

	/**
	 * Transforms the given vector of complex values into their inverse Fourier
	 * Transform. The resulting vector is a vector of doubles only (the real
	 * part of the inverse transformation). The vector of complex values must
	 * be the hermitian redundant version of the FFT transform.
	 *
	 * @param input The vector of complex numbers to transform
	 * @param output The container of doubles where the result will be stored.
	 */
	template <typename T>
	void backward(const std::vector<std::complex<double>> &input, T &output) const;

	unsigned int get_size() {
		return size;
	}

	unsigned int get_hermitian_size() {
		return hermitian_size;
	}

private:
	unsigned int size;
	unsigned int hermitian_size;
	unsigned int omp_threads;
	std::unique_ptr<double, fftw_deleter<double>> real_buf;
	std::unique_ptr<fftw_complex, fftw_deleter<fftw_complex>> complex_buf;
	fftw_plan forward_plan;
	fftw_plan backward_plan;
};

}  // namespace profit

#endif /* PROFIT_FFWT */

#endif // PROFIT_FFT_IMPL_H