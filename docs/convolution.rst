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

Supported convolution methods
-----------------------------

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
* :class:`OpenCLConvolver` is a brute-force convolver
  implemented in OpenCL.
  It offers both single and double floating-point precision
  and its performance is usually better
  that that of the `BruteForceConvolver`.

Creating a Convolver
--------------------

Instead of manually selecting the class that should be used,
users create :class:`Convolver` instances
via the :func:`create_convolver` function.
:func:`create_convolver` lets the user specify
which type of convolver should be created
(either using an enumeration, or a standard string value),
and a set of creation preferences
that apply differently to different types of Convolvers.
Once created,
users can call the :member:`Convolver::convolve` method
directly on the resulting convolver,
or assign it to a :class:`Model` instance
for it to use it.

If a :class:`Model` needs to perform convolution
and a :class:`Convolver` has been set
on its :member:`Model::convolver` member
then that convolver is used.
If no convolver has been set,
it creates a new :class:`BruteForceConvolver`
and uses that to perform the convolution.
