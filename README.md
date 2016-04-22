# libprofit

libprofit is a C library that produces images based on different luminosity profiles.

The currently supported profiles (and their parameters) are:

 * **sersic**: X/Y center, sersic number, magnitude, angle, effective radius, axes ratio, boxing
 * **sky**: background

With time we intend to add more profiles to the library. Users can also provide their own profiles.

Any number of profiles can be given, all of which are calculated and summed up, resulting in a final image.

libprofit has no compiling dependencies other than libc and libm; however some profiles require some high-level functions that can be found in third-party mathematical packages like R or GSL. These functions must thus be provided when linking libprofile into the resulting library or program.
