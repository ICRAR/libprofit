Using libprofit
===============

.. default-domain:: c

At the core of libprofit sits :type:`profit_model`.
This structure holds all the information needed to generate an image.
Different profiles (instances of :type:`profit_profile`)
are appended to the model, which is then evaluated.

The basic usage pattern then is as follows:

#. First obtain a model instance::

	 profit_model *model = profit_create_model();

#. Create a profile (for a list of supported names see :doc:`profiles`)::

	 profit_profile *sersic_profile = profit_create_profile("sersic");

#. Customize your profile.
   An explicit cast must be performed on the :type:`profit_profile` to turn it
   into the specific profile type.
   By convention these sub-types are named after the profile they represent,
   like this::

	 profit_sersic_profile *sp = (profit_sersic_profile *)sersic_profile;
	 sp->xcen = 34.67;
	 sp->ycen = 9.23;
	 sp->axrat = 0.345;
	 [...]

#. Add the new profile to the model::

	 profit_add_profile(model, sersic_profile);

#. Repeat the previous three steps for all profiles
   you want to include in your model.

#. Evaluate the model simply run::

	 profit_eval_model(model);

#. After running check if there are have been errors
   while generating the image.
   If no errors occurred you can safely access the data
   stored in :member:`profit_model.image`::

	 char *error = profit_get_error(model);
	 if( error ) {
	     printf("Oops! There was an error evaluating the model: %s", error);
	 }
	 else {
	    do_something_with_your_image(model->image);
	 }

#. Finally dispose of the model.
   This should **always** be called,
   regardless of whether the model was actually used or not,
   or whether its evaluation was successful or not::

	 profit_cleanup(model);

