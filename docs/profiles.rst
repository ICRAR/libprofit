Profiles
========

This is the list of profiles supported by libprofit.

.. contents::

Sersic
------

The sersic profile follows I don't know what.

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
