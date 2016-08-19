Getting *libprofit*
###################

*libprofit* is currently hosted in `GitHub <https://github.com/rtobar/libprofit>`_.
To get a copy you can clone the repository::

 git clone https://github.com/rtobar/libprofit


Compiling
=========

*libprofit* depends (optionally) on:

* `GSL <https://www.gnu.org/software/gsl/>`_
* `R <https://www.r-project.org/>`_

Both dependencies satisfy the same requirements,
so they are mutually exclusive.
If both are present GSL takes precedence.
If none is present,
the ``sersic`` profile will not work.

*libprofit*'s compilation system is based
on `cmake <https://cmake.org/>`_.
``cmake`` will check that you have a proper compiler
(anything supporting some basic C++11 should do),
and scan the system for dependencies.

To compile *libprofit* run
(assuming you are inside the ``libprofit`` directory already)::

 $> mkdir build
 $> cd build
 $> cmake ..
 $> make
 $> # optionally for system-wide installation: sudo make install

With ``cmake`` you can also specify additional compilation flags.
If you want to generate the fastest possible code try this::

 $> CXXFLAGS="-O3 -march=native" cmake ..

You can also specify a different installation directory like this::

 $> cmake -DCMAKE_INSTALL_PREFIX=~/my/installation/directory

Please refer to the ``cmake`` for further options.
