Changelog
=========

.. default-domain:: cpp
.. highlight:: cpp
.. namespace:: profit

.. rubric:: Development version

* :class:`Model` offers a new flavour of the ``evaluate`` function
  where users can give a pre-existing :class:`Image` object
  on which the result will be written.
  If this user-provided :class:`Image` object
  is not of the size internally needed by :class:`Model`
  then it will be internally resized automatically.
* A new ``get_drawing_dimensions`` method
  has been added to the :class:`Model` class.
* A new ``null`` convolver has been added
  that does no convolution
  and simply returns the source image unmodified.
  This is only useful for testing.

.. rubric:: 1.9.3

* A bug in the OpenCL implementation of the radial profiles
  prevented Models with multiple profiles
  from displaying correctly,
  as the output image would contain
  only the values of last profile.
  This was a problem introduced
  only in the last version of *libprofit*,
  and not an ongoing issue.
* When using OpenCL,
  any radial profile specifying ``rough=true``
  caused the output image not to be scaled properly,
  with values not taking into account the profile's magnitude
  or pixel scale.
  This seems to have been an issue for a long time,
  but since ``rough=true`` is not a common option
  it had gone under the radar for some time.

.. rubric:: 1.9.2

* All profile evaluation has been changed
  from being absolute (profiles set the final value of a pixel)
  to be additive (their add their pixel values onto the image).
  This change in behavior has the effect
  that one less memory allocation is needed,
  which can be a big difference
  when generating big images,
  while also simplifying the logic
  of the :class:`Model` evaluation.
* :class:`Model` objects now internally store
  the normalized version of the PSF
  given by the user instead of the original,
  which was never really needed.
* :program:`profit-cli` now makes it easier to specify
  multiple copies of the same profile,
  useful for scaling tests.
  Also, writing FITS files in little endian systems
  doesn't allocate extra memory anymore.
* Minor improvements to imaging classes.

.. rubric:: 1.9.1

* The implementation of the :class:`Model` class has been improved.
  In particular it has been made more memory efficient,
  which is particularly important in scenarios
  where many profiles (in the order of thousands)
  are added into it.
  Previously each profile was allocated its own :class:`Image`,
  which added both to the memory footprint,
  and to the total runtime.
  Now a single scratch space is used for all profiles,
  and individual results are immediately summed up,
  respecting the convolution settings of each profile.
  Experiments with the :ref:`null profile <profiles.null>`
  show a significant decrease in runtime
  when many Model evaluations take place.

.. rubric:: 1.9.0

* Implemented correct :doc:`flux capturing <flux_capturing>`.
  This feature was previously implemented
  in the `ProFit <https://github.com/ICRAR/ProFit>`_ R package
  as part of its fitting process,
  but it was otherwise unavailable.
* Added explicit support to allow convolution
  of images against kernels with bigger dimensions
  than the images themselves.
  This was previously supported implicitly, and only in certain cases,
  by the OpenCL convolver,
  while the FFT convolver threw an proper exception,
  and the brute-force convolvers usually crashed.
  This first implementation is not ideal,
  but the use case is rare.
* Several performance and code improvements,
  like removing unnecessary code,
  avoiding unnecessary conversions
  and avoiding a few dynamic allocations.

.. rubric:: 1.8.2

* Users can now select the underlying
  SIMD-capable instruction set to use
  for brute-force convolution.
* New library method :func:`has_simd_instruction_set`
  for users to check whether libprofit was compiled
  with support for different instruction sets.
* Improved FFTW-based convolver performance
  by avoiding dynamic memory allocation at convolution time.
  This brings a noticeable performance improvement
  of around 20%.

.. rubric:: 1.8.1

* Adding support for FFTW versions lower than 3.3.

.. rubric:: 1.8.0

* :program:`profit-cli` compiling in Windows.
* New :func:`Profile::parameter` method to specify
  parameters and values with a single ``name=value`` string.
* New utility methods: :func:`trim`, :func:`split` and :func:`setenv`.
* Using SSE2/AVX SIMD extensions to implement brute-force convolution
  if the CPU supports it, with pure C++ implementation as a fallback.
  Can be disabled with ``-DLIBPROFIT_NO_SIMD=ON``.
* Potentially fixed the importing of FFTW wisdom files
  in systems with more than one FFTW installation.
* Fixed compilation of ``brokenexponential`` OpenCL kernel in platforms where it
  was failing to compile.
* Compiling in release mode (i.e., ``-O3 -DNDEBUG`` in gcc/clang) by default.
* Lowering OpenMP requirement to 2.0 (was 3.0).
* OpenCL kernel cache working for some platforms/devices that was not
  previously working.
* Many internal code cleanups and design changes
  to make code easier to read and maintain.

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

* Added :func:`Image::upsample` and :func:`Image::downsample` to scale an
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
