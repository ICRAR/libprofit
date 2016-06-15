Getting *libprofit*
###################

*libprofit* is currently hosted in `GitHub <https://github.com/rtobar/libprofit>`_.
To get a copy you can clone the repository::

 git clone https://github.com/rtobar/libprofit

*libprofit*'s compilation system is still extremely simple.
It uses ``make`` only to drive the compilation process,
leaving the user with the responsibility of adjusting
the ``CC``, ``CFLAGS`` and ``LDFLAGS`` variables as needed.
Compiling it then is as simple as::

 $> cd src/
 $> make

*libprofit* currently understands the following pre-processing macros:

* ``HAVE_GSL``: Indicates that the GSL is available in the system.
  If given then GSL's implementation of a handful of statistic functions
  like the beta and gamma functions is used in the code.
  Note that if you pass down this option
  then the corresponding linking options must be given as well
  to link *libprofit* with the GSL library.
* ``HAVE_R``: Indicates that the R library is available in the system.
  If given then R's implementation of the functions mentioned avobe
  is used instead.
  Again, passing down this option means
  that the corresponding linking option must be given as well.

For example, if compiling *libprofit* with built-in support for the GSL
then the following command will do::

 CFLAGS=-DHAVE_GSL LDFLAGS="-lgls -lgslcblas" make

In the future we might support a more automatic compilation system
like ``cmake`` or ``autotools``.
