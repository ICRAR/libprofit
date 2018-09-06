Changelog
=========

.. default-domain:: cpp
.. highlight:: cpp
.. namespace:: profit

.. rubric:: Development version

* :program:`profit-cli` compiling in Windows.
* New :func:`Profile::parameter` method to specify
  parameters and values with a single ``name=value`` string.
* New utility methods: :func:`trim`, :func:`split` and :func:`setenv`.
* Fixed compilation of ``brokenexponential`` OpenCL kernel in platforms where it
  was failing to compile
* Lowering OpenMP requirement to 2.0 (was 3.0)
* OpenCL kernel cache working for some platforms/devices that was not
  previously working

.. rubric:: 1.7.4

* FFT convolution using hermitian redundancy. This increases performance of
  FFT-based convolution by at least 10% in release builds, and addresses some
  warnings pointed out by ``valgrind``.

.. rubric:: 1.7.3

* Added :func:`init_diagnose` and :func:`finish_diagnose` functions to avoid
  printing to stdout/stderr from within libprofit.

.. rubric:: 1.7.2

* Fixed ``double`` detection support for OpenCL devices regardless of the
  supported OpenCL version.
* Fixed a few compiling issues under Visual Studio compiler.
* Continuous integration in Windows via `AppVeyor <https://ci.appveyor.com/project/rtobar/libprofit>`_

.. rubric:: 1.7.1

* Added :member:`Image::upsample` and :member:`Image::downsample` to scale an
  image up or down (using different modes).
* Added :func:`Model::set_return_finesampled` to return internally
  upsampled images.

.. rubric:: 1.7.0

* Internal implementation dependencies clearly hidden from users. This means
  that users compiling against libprofit don't need to search for header files
  other than libprofit's, making it much easier to write code against libprofit.
* :class:`Model` redesigned. No member variables are exposed anymore; instead
  different setter/getter methods must be used.
* :class:`Image` redesigned. In summary, it looks much more like a standard
  container now.
* New :func:`Model::set_crop` specifies whether cropping should be carried out
  after convolution, if the convolution needs to pad the image.
* :func:`Model::evaluate` has an extra optional parameter to receive the
  offset at which cropping needs to happen (if it hasn't, see
  :func:`Model::set_crop`) to remove padding from the resulting image.
* FFTW convolution uses real-to-complex and complex-to-real forward and
  backwards transforms respectively (instead of complex-to-complex transforms
  both ways), which is more efficient and should use less memory.
* New on-disk OpenCL kernel cache. This speeds up the creation of OpenCL
  environments by a big factor as compilation of kernels doesn't happen every
  time an environment is created.
* New on-disk FFTW plan cache. This speeds up the creation of FFT-based
  convolvers by a big factor as the plans are not calculated every time for a
  given set of parameters.
* New ``null`` profile, useful for testing.
* New :func:`init` and :func:`finish` calls to initialize and finalize
  libprofit. These are mandatory, and should be called before and after using
  anything else from libprofit.

.. rubric:: 1.6.1

* Brute-force convolver about 3x faster than old version.
* Fixing compilation failure on MacOS introduced in 1.6.0.
* Center pixel in sersic profile treated specially only if ``adjust`` parameter
  is on.
