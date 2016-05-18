API
===

.. default-domain:: c

.. type:: profit_model

   The root structure holding all the information needed by libprofit
   to generate an image.
   Users should create profit_model instances via :func:`profit_create_model`
   and destroy them using :func:`profit_cleanup`.

.. member:: unsigned int profit_model.width

   The width of the image that profit will generate for this model.
   It must be greater than 0.

.. member:: unsigned int profit_model.height

   The height of the image that profit will generate for this model.
   It must be greater than 0.

.. function:: profit_model* profit_create_model(void)

   Returns a newly created model to be used for evaluation.

.. function:: profit_profile* profit_create_profile(const char *profile_name)

   Returns a newly created profile that corresponds to profile_name.

.. function:: void profit_cleanup(profit_model *model)

   Releases all resources associated with the given model, after which
   it cannot be used anymore.
