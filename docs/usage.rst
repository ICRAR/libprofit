Using *libprofit*
=================

From the command-line
---------------------

*libprofit* ships with a command-line utility ``profit-cli``
that reads the model and profile parameters from the command-line
and generates the corresponding image.
It supports all the profiles supported by *libprofit*,
and can output the resulting image as text values, a binary stream,
or as a simple FITS file.

Run ``profit-cli -h`` for a full description on how to use it,
how to specify profiles and model parameters,
and how to control its output.

Programatically
---------------

As it name implies, *libprofit* also ships a shared library
exposing an API that can be used by any third-party application.
This section gives a brief overview on how to use this API.
For a full reference please refer to :doc:`api`.

.. default-domain:: cpp
.. highlight:: cpp
.. namespace:: profit

At the core of *libprofit* sits :class:`Model`.
This class holds all the information needed to generate an image.
Different profiles (instances of :class:`Profile`)
are appended to the model, which is then evaluated.

The basic usage pattern then is as follows:

#. Add the profit include::

	 #include <profit/profit.h>

#. Initialize the library with the :func:`init` function.
   This needs to be called *only once* in your program::

	 profit::init();

#. First obtain a model instance that will generate profile images
   for a given width and height::

	 profit::Model model(width, height);

#. Create a profile. For a list of supported names see :doc:`profiles`;
   if you want to support a new profile see :doc:`new_profile`.
   If an unknown name is given an :class:`invalid_parameter` exception will be
   thrown::

	 profit::ProfilePtr sersic_profile = model.add_profile("sersic");

#. Customize your profile.
   To set the different parameters on your profile call
   :func:`Profile::parameter` with the parameter name and value::

	 sersic_profile.parameter("xcen", 34.67);
	 sersic_profile.parameter("ycen", 9.23);
	 sersic_profile.parameter("axrat", 0.345);
	 sersic_profile.parameter("nser=3.56");
	 // ...

   A complete list of parameters can be found on and :doc:`profiles` and
   :doc:`api`.

#. Repeat the previous two steps for all profiles
   you want to include in your model.

#. Evaluate the model simply run::

	 profit::Image result = model.evaluate();

#. If the resulting image needs to be cropped
   (see :ref:`convolution.image_cropping` for full details)
   an additional argument needs to be passed
   to :func:`Model::evaluate`
   to receive the offset at which cropping needs to be,
   like this::

	 profit::Point offset;
	 profit::Image result = model.evaluate(offset);
	 profit::Image cropped_image = result.crop({width, height}, offset);

#. If there are have been errors
   while generating the image
   an :class:`invalid_parameter` exception will be thrown by the code,
   so users might want to use a ``try/catch`` statement
   to identify these situations::

	 try {
	     auto result = model.evaluate();
	 } catch (profit::invalid_parameter &e) {
	     cerr << "Oops! There was an error evaluating the model: " << e.what() << endl;
	 }

#. When the model is destroyed the underlying profiles are destroyed as well.

#. When you are finished using the library,
   call the :func:`finish` function::

	 profit::finish();

To illustrate this process, refer to the following figure:

.. image:: images/evaluation.png
