API
===

.. contents::

.. default-domain:: cpp
.. namespace:: profit

*libprofit* has two main data types:
the model class (:class:`Model`)
and the base profile class (:class:`Profile`).
We introduce the base profile first, then the model.

As a small reference, see the following class diagram.

.. image:: images/types.png

Exceptions
----------

.. class:: invalid_parameter

   An exception thrown by :class:`Profile` and :class:`Model` objects
   to signal that an invalid parameter was set on them,
   preventing an image from being calculated.

Profile class
-------------

.. class:: Profile

   The base class that ever profile class must inherit from.
   It contains all the shared aspects across profiles,
   like a name. It also specifies the methods
   that validate and evaluate a profile, and that should be
   implemented by each profile subclass.

   Users should create :class:`Profile` instances
   via :member:`Model::add_profile`.

.. namespace-push:: Profile

.. function:: virtual void validate() = 0

   Abstract method to be implemented by subclasses.
   It checks that the parameters supplied to the profile are valid,
   and throws an :class:`invalid_parameter` exception otherwise.

.. function:: virtual void evaluate(double * image) = 0

   Abstract method to be implemented by subclasses.
   It evaluates the profile and stores the result in ``image``.

.. member:: std::string name

   The name of this profile.

.. member:: bool convolve

   A boolean flag indicating whether the image produced by this profile
   should be convolved with the model's PSF or not.
   Setting this flag to ``true`` but failing to provide a PSF
   results in an error.

.. namespace-pop::

Model class
-----------

.. class:: Model

   The root object holding all the information needed by *libprofit*
   to generate an image.

.. namespace-push:: Model

.. function:: Profile * add_profile(std::string profile_name)

   Creates a new profile for the given name and adds it to the given model.
   On success, the new profile is created, added to the model,
   and its reference is returned for further customization.
   On failure (i.e., if a profile with the given name is not supported)
   ``NULL`` is returned and no profile is added to the model.

.. function:: void evaluate()

   Calculates an image using the information contained in the model.
   The result of the computation is stored in the image field.

.. member:: unsigned int width

   The width, in pixels, of the image that profit will generate for this model.
   It must be greater than 0.
   See :doc:`coordinates` for more details.

.. member:: unsigned int height

   The height, in pixels, of the image that profit will generate for this model.
   It must be greater than 0.
   See :doc:`coordinates` for more details.

.. member:: double * image

   The image produced by this model.
   The image has the dimensions specified in the model.
   If there was any error when evaluating the model
   this field will remain unset.

.. member:: unsigned int res_x

   The span of the horizontal coordinate of the image that profit will generate
   for this model.
   It must be greater than 0.
   See :doc:`coordinates` for more details.

.. member:: unsigned int res_y

   The span of the vertical image coordinate.
   It must be greater than 0.
   See :doc:`coordinates` for more details.

.. member:: double magzero

   The zero magnitude of this model.

.. member:: std::vector<Profile *> profiles

   A vector of pointers to the individual profiles
   used to generate the model's image.

.. member:: double * psf

   An array containing the values of a Point Spread Function (PSF).
   The PSF is used to convolve the profiles that request convolving,
   and as the source image of the ``psf`` profile.

.. member:: unsigned int psf_width

   The width of the PSF image.

.. member:: unsigned int psf_height

   The height of the PSF image.

.. member:: bool * calcmask

   A boolean mask with the same dimensions of the model
   that indicates for each pixel of the image
   whether the profiles should be calculated or not.
   If ``NULL`` all pixels are calculated.
