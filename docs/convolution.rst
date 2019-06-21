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
(via their :func:`Convolver::convolve` method).
Depending on the size of the problem,
and on the libraries available on the system,
different convolver types will be available to be used:

* :enumerator:`BRUTE_OLD` is the simplest convolver.
  It implements a simple, brute-force 2D convolution algorithm.
* :enumerator:`BRUTE` is a brute-force convolver
  that performs better that :enumerator:`BRUTE_OLD`, but still
  implements simple, brute-force 2D convolution. It is the default
  convolver used by a :class:`Model` that hasn't been assigned one,
  but requires one.
* :enumerator:`FFT` is a convolver
  that uses Fast Fourier transformations to perform convolution.
  Its complexity is lower than the :enumerator:`BRUTE`,
  but its creation can be more expensive.
* :enumerator:`OPENCL` is a brute-force convolver
  implemented in OpenCL.
  It offers both single and double floating-point precision
  and its performance is usually better
  that that of the :enumerator:`BRUTE`.

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

If a :class:`Model` needs to perform convolution
and a :class:`Convolver` has been set
on its :member:`Model::convolver` member
then that convolver is used.
If no convolver has been set,
it creates a new :enumerator:`BRUTE`
and uses that to perform the convolution.


Using a convolver
-----------------

Once created,
users can call the :func:`Convolver::convolve` method
directly on the resulting convolver,
(or assign it to a :class:`Model` instance for it to use it).
The :func:`Convolver::convolve` methods needs at least three parameters:
an image, a kernel and a mask.
Convolvers will convolve the image with the kernel
only for the pixels in which the mask is set,
or for all pixels if an empty mask is passed.
This implies that the mask, if not empty,
must have the same dimensions that the image.


.. _convolution.image_cropping:

Image cropping
--------------

Some convolvers internally work
with images that are larger
than the original source image
(mostly due to efficiency reasons).
After this internal image expansion occurs,
and the convolution takes place,
the resulting image
is usually cropped at the corresponding point
to match original source image size and positioning
before being returned to the user.

However, users might want to pick
into this internal, non-cropped result
of the convolution process.
To do this,
an additional ``crop`` parameter
in the :func:`Convolver::convolve` method
determines whether the convolver should return
the original, and potentially bigger, image.
When a non-cropped image is returned,
an additional ``offset_out`` parameter
can be given to find out the offset
at which cropping would have started.
The cropping dimensions do not need to be queried,
as they always are the same
of the original source image given to the convolver.


.. _convolution.model:

Model convolution
-----------------

During model evaluation (i.e., a call to :func:`Model::evaluate`)
users might want to be able to retrieve the non-cropped result
of the internal convolution that takes place
during model evaluation
(as explained in :ref:`convolution.image_cropping`).

To do this, users must first
call :func:`Model::set_crop` with a ``false`` argument.
When calling :func:`Model::evaluate`,
users must then also give a :type:`Point` argument
to retrieve the offset at which
cropping should be done
to remove the image padding
added by the convolution process.
