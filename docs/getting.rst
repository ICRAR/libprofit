Getting *libprofit*
###################

*libprofit* is currently hosted in `GitHub <https://github.com/ICRAR/libprofit>`_.
To get a copy you can clone the repository::

 git clone https://github.com/ICRAR/libprofit


Compiling
=========

*libprofit* depends on:

* `GSL <https://www.gnu.org/software/gsl/>`_
* `R <https://www.r-project.org/>`_

Both dependencies satisfy the same requirements,
so they are mutually exclusive,
but at least one of them is necessary.
If both are present GSL takes precedence.

Optional requirements are:

* An `OpenMP <http://www.openmp.org/>`_-enabled compiler
* An `OpenCL <https://www.khronos.org/opencl/>`_ installation
* `FFTW <http://www.fftw.org/>`_

*libprofit*'s compilation system is based
on `cmake <https://cmake.org/>`_.
``cmake`` will check that you have a proper compiler
(anything supporting some basic C++11 should do),
and scan the system for all required dependencies.

To compile *libprofit* run
(assuming you are inside the ``libprofit`` directory already)::

 $> mkdir build
 $> cd build
 $> cmake ..
 $> make
 $> # optionally for system-wide installation: sudo make install

With ``cmake`` you can also specify additional compilation flags.
For example, if you want to generate the fastest possible code
you can try this::

 $> cmake .. -DCMAKE_CXX_FLAGS="-O3 -march=native"

You can also specify a different installation directory like this::

 $> cmake .. -DCMAKE_INSTALL_PREFIX=~/my/installation/directory

Other ``cmake`` options that can be given in the command-line include:

* ``LIBPROFIT_USE_R``: prefer R libraries over GSL libraries
* ``LIBPROFIT_TEST``: enable compilation of unit tests
* ``LIBPROFIT_DEBUG``: enable debugging-related code
* ``LIBPROFIT_NO_OPENCL``: disable OpenCL support
* ``LIBPROFIT_NO_OPENMP``: disable OpenMP support
* ``LIBPROFIT_NO_FFTW``: disable FFTW support

Please refer to the ``cmake`` documentation for further options.
