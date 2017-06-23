Convolution
===========

.. contents:: Contents
   :local:

.. default-domain:: cpp
.. highlight:: cpp
.. namespace:: profit

Image convolution in *libprofit* happens optionally
as part of a :class:`Model` evaluation.
Internally, the :class:`Model` uses a :class:`Convolver`
to perform convolution.

Convolvers are objects that carry out convolution
(via their :member:`Convolver::convolve` method).
Depending on the size of the problem,
and on the libraries available on the system,
different convolver types will be available to be used:

* :class:`BruteForceConvolver` is the simplest convolver.
  It implements a simple, brute-force 2D convolution algorithm.
* :class:`FFTConvolver` is a convolver
  that uses Fast Fourier transformations to perform convolution.
  Its complexity is lower than the :class:`BruteForceConvolver`,
  but its creation can be more expensive.

If a :class:`Model` needs to perform convolution
and a :class:`Convolver` has been set
on its :member:`Model::convolver` member
then that convolver is used.
If no convolver has been set,
it creates a new one
depending on its settings
(namely, on the value of
:member:`Model::use_fft`, :member:`Model::omp_threads`,
:member:`Model::fft_effort` and :member:`Model::reuse_psf_fft`),
and stores it for later retrieval.
