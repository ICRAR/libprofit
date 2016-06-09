API
===

.. default-domain:: c

Data types
----------

*libprofit* has two main data types:
the model structure (:type:`profit_model`)
and the base profile structure (:type:`profit_profile`).
We introduce the base profile first, then the model.

.. type:: profit_profile

   The base structure that ever profile type must use.
   It contains all the shared aspects across profiles,
   like a name, an error string, and pointers to the functions
   that initialize and evaluate a profile

   A member of this type **must** be declared as the first structure member
   of each profile type structure.
   This allows them to share the same memory address,
   and therefore make a pointer safely castable from one type to the other.

   Users should create :type:`profit_profile` instances
   via :func:`profit_create_profile`,
   and add them to the corresponding model via :func:`profit_add_profile`.

.. member:: char * profit_profile.name

   The name of this profile.

.. member:: char * profit_profile.error

   An error string indicating that an error related to this profile was
   detected. The error string can be set either during the profile
   initialization or during the image creation process. Users should check
   that there is no error in any of the profiles after making a model
   using :func:`profit_get_error`.

.. member:: bool profit_profile.convolve:

   A boolean flag indicating whether the image produced by this profile
   should be convolved with the model's PSF or not.
   Setting this flag to ``true`` but failing to provide a PSF
   results in an error.



.. type:: profit_model

   The root structure holding all the information needed by *libprofit*
   to generate an image.
   Users should create :type:`profit_model` instances via :func:`profit_create_model`
   and destroy them using :func:`profit_cleanup`.

.. member:: unsigned int profit_model.width

   The width, in pixels, of the image that profit will generate for this model.
   It must be greater than 0.
   See :doc:`coordinates` for more details.

.. member:: unsigned int profit_model.height

   The height, in pixels, of the image that profit will generate for this model.
   It must be greater than 0.
   See :doc:`coordinates` for more details.

.. member:: double * profit_model.image

   The image produced by this model.
   The image has the dimensions specified in the model.
   Users should check if there was any error when evaluating the model
   using :func:`profit_get_error`, in which case this field will remain unset.

.. member:: unsigned int profit_model.res_x

   The span of the horizontal coordinate of the image that profit will generate
   for this model.
   It must be greater than 0.
   See :doc:`coordinates` for more details.

.. member:: unsigned int profit_model.res_y

   The span of the vertical image coordinate.
   It must be greater than 0.
   See :doc:`coordinates` for more details.

.. member:: double profit_model.magzero

   The zero magnitude of this model.

.. member:: unsigned int profit_model.n_profiles

   The number of profiles used to generate the model's image.

.. member:: profit_profile ** profit_model.profiles

   A list of pointers to the individual profiles
   used to generate the model's image.

.. member:: double * profit_model.psf

   An array containing the values of a Point Spread Function (PSF).
   The PSF is used to convolve the profiles that request convolving,
   and as the source image of the ``psf`` profile.

.. member:: unsigned int profit_model.psf_width

   The width of the PSF image.

.. member:: unsigned int profit_model.psf_height

   The height of the PSF image.

.. member:: bool * profit_model.calcmask

   A boolean mask with the same dimensions of the model
   that indicates for each pixel of the image
   whether the profiles should be calculated or not.
   If ``NULL`` all pixels are calculated.

.. member:: char * profit_model.error

   An error string indicating that an error at the model level has been
   detected.
   Users should check that there is no error in any of the profiles
   after making a model using :func:`profit_get_error`.

Functions
---------

These are the set of functions
that are externally visible from *libprofit* to the users.
For an example on how to use them see :doc:`usage`.

.. function:: profit_model * profit_create_model(void)

   Creates a new model to which profiles can be added, and that can
   be used to calculate an image.

.. function:: profit_profile * profit_create_profile(const char * profile_name)

   Creates a new profile for the given name.
   On success, the new profile is created and its reference is returned for
   further customization.
   On failure (i.e., if a profile with the given name is not supported)
   ``NULL`` is returned.

.. function:: void profit_add_profile(profit_model * model, profit_profile * profile)

   Adds the given profile to the model.

.. function:: void profit_eval_model(profit_model * model)

   Calculates an image using the information contained in the model.
   The result of the computation is stored in the image field.

.. function:: char * profit_get_error(profit_model * model)

   Returns the first error string found either on the model itself or in any of
   it profiles. This method should be called on the model right after invoking
   profit_eval_model to make sure that no errors were found during the process.
   If ``NULL`` is returned it means that no errors were found and that the image
   stored in the model is valid.

.. function:: void profit_cleanup(profit_model * model)

   Frees all the resources used by the given model, after which it cannot be
   used anymore.
