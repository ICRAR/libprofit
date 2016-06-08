Profiles
========

.. contents::

This section lists the profiles currently supported by *libprofit*.

``sersic``
----------

An implementation of the
`Sersic luminosity profile <https://en.wikipedia.org/wiki/Sersic_profile>`_.
The sersic profile describes the intensity of a galaxy
depending on its distance to the center.

The sersic profile accepts the following parameters:

* **xcen**: x centre of the Sersic profile (can be fractional pixel positions).
* **ycen**: y centre of the Sersic profile (can be fractional pixel positions).
* **mag**: Total magnitude of the Sersic profile.
  Converted to flux using ``flux = 10^(-0.4*(mag - magzero))``, where
  ``magzero`` is that of the containing model.
* **re**: Effective ratio
* **nser**: Sersic index of the Sersic profile.
* **ang**: The orientation of the major axis of the Sersic profile, in degrees.
  The starting point is the positive Y image axis and grows counter-clockwise.
* **axrat**: Axial ratio of the Sersic profile defined as minor-axis/major-axis,
  i.e. ``axrat = 1`` is a circle and ``axrat = 0`` is a line.
* **box**: The boxiness of the Sersic profile that traces contours of iso-flux,
  defined such that r = (x^(2+box)+y^(2+box))^(1/(2+box)).
  When ``box = 0`` the iso-flux contours will be normal ellipses,
  but modifications between ``-1 < box < 1`` will produce visually boxy distortions.
  Negative values have a pin-cushion effect, whereas positive values have a barrel effect
  (the major and minor axes staying fixed in all cases).

The sersic profile implements recursive sub-pixel sampling for better results
in areas closer to the profile center.
This sub-sampling can be controller by the following additional parameters:

* **rough**: Don't perform any sub-sampling, ignore the rest of the parameters.
* **re_switch**: Effective radius fraction within which sub-sampling is performed.
  Pixels outside this radius are not sub-sampled.
* **max_recursions**: The maximum levels of recursions allowed.
* **resolution**: Resolution (both horizontal and vertical) to be used
  on each new recursion level.
* **acc**: Accuracy after which recursion stops.

The sersic profile also implements far-pixel filtering,
quickly zeroing pixels that are too far away
from the profile center.
This filtering can be controller by the following parameters:

* **re_max**: Maximum *re*-based distance to consider for filtering.
* **rescale_flux**: Whether the calculated profile flux should be scaled
  to take into account the filtering performed by **re_max**.
* **calcmask**: A logical 2D map that indicates whether the profile
  should be calculated at a given pixel or not.
  By default all pixels are calculated.

Finally, an **adjust** parameter allows the user
whether adjustments of most of the parameters described
above should be done automatically depending on the profile parameters.
*libprofit* makes a reasonable compromise between speed and accuracy,
and therefore this option is turned on by default.

``sky``
-------

The sky profile provides a constant value for an entire image.

* **bg**: Value per pixel for the background.
  This should be the value as measured in the original image,
  i.e. there is no need to worry about the effect of the model's ``magzero``.

``psf``
-------

The psf profile adds the model's psf to the model's image
at a specific location and for a given user-defined magnitude.

* **xcen**: The x position at which to generate the centred PSF
  (can be fractional pixels).
* **ycen**: The y position at which to generate the centred PSF
  (can be fractional pixels).
* **mag**: The total flux magnitude of the PSF.
