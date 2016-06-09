Adding a profile
================

.. default-domain:: c
.. highlight:: c

.. contents::

This section explains
the steps required to add a new profile to *libprofit*.

A profile needs three parts to be complete:

* A new data type that contains
  all the information needed by the profile
* A set of functions that create it and evaluate it
* A association with a standard name

In all steps below,
a completely artificial ``example`` profile is being added,
This new profile takes three parameters:
``param1`` and ``param2`` are double numbers,
while ``param3`` is an integer.
The profile fills the image by taking the X and Y coordinates
and filling the pixel
with the value ``|(param1 - param2) * param3 * (x - y)|``
and requires that all parameters are positive or 0.

The data types used in this example
are described in detail in :doc:`api`.

Data type
---------

The first step to add a new profile is to
define the C structure that will hold all its information.
Any kind of information can be added to the structure,
but it is **required** that the first member of the structure
is of type ``profit_profile``.
This member is used internally by *libprofit*
to handle all profile types as being of the same kind.
The current profiles name this first member simply ``profile``,
so the same name is suggested for new profiles.
It is suggested you also ``typedef`` your new structure for easier
reference in the future.
Again, any name can be chosen,
but to keep things consistent
``profit_xxxx_profile`` is used throughout the code,
with ``xxxx``` being the name of the profile.

The structure above should be defined
in an ``.h`` file in the ``include`` directory
so it can be used by others.

All in all, it should look like this::

 typedef struct _profit_example_profile {
     profit_profile profile;
     double param1;
     double param2;
     int param3;
 } profit_example_profile;

Functions
---------

Each profile requires a minimum of three functions
that need to be written

* One to create a new profile,
* One to initialize it, and
* One to evaluate it

They can carry any name,
but the current convention is that they are called
``profit_create_xxxx``, ``profit_init_xxxx`` and ``profit_make_xxxx``
respectively, where ``xxxx`` is the name of the profile.
The functions should live all
in the same file ``.c`` so they can easily see each other.

Initialization
^^^^^^^^^^^^^^

We will start by looking at the initialization function.
Its signature looks like this::

 void profit_init_example(profit_profile *profile, profit_model *model);

The initialization function takes in two arguments:
the containing ``model``
and the particular ``profile`` being initialized
(see :doc:`structure`).
Its responsibility is to validate the inputs of the profile.

The first thing you want to do in the function
is to cast the profile into the ``profit_example_profilt`` type
to have access to its specific members.
In the case of the ``example`` profile it was mentioned
that all parameters must be positive,
so the code must test for that.
If a violation occurs, an error string is recorded
n the ``profile->error`` field.
This error will prevent the profile (and in fact the whole model)
from being evaluated.

An example implementation would thus look like this::

 void profit_init_example(profit_profile *profile, profit_model *model) {

     profit_example_profile *e = (profit_example_profile *)profile;

     if ( e->param1 < 0 ) {
        profile->error = strdup("param1 is negative");
        return;
     }
     if ( e->param1 < 0 ) {
        profile->error = strdup("param2 is negative");
        return;
     }
     if ( e->param3 < 0 ) {
        profile->error = strdup("param3 is negative");
        return;
     }

 }

Note that having access to the model means
that one can validate profile-specific values
against model-global values as well.
For example, if a new restriction is added stating
that the ``example`` profile can only be run on images
that are bigger than 20 x 20
then the following code could be added::

 if ( model->width < 20 || model->height < 20 ) {
     profile->error = strdup("can't apply example profile to images less than 20x20");
     return;
 }

Finally, if a profile needs no validation at all
an initialization function must still be provided
with an empty body.

Evaluation
^^^^^^^^^^

Next, we loop to the evaluation profile function.
Its signature looks like this::

 void profit_make_example(profit_profile *profile, profit_model *model, double *image);

Just like in the initialization function,
both the containing model
and the specific profile to evaluate
are given,
and probably the first thing you want to do
is to cast the profile into the ``profit_example_profilt`` type
to have access to its specific members.

A third ``image`` argument is also received.
This corresponds to the surface where the pixels must be drawn.
Each profile in the model receives a different image surface,
so you will start with a clean slate.
The image is already initialized with zeros,
so if your profile doesn't cover the entire image
there is no need to cover unfilled areas with zeros.
The image is an array of ``double`` values
representing a 2D image with dimensions ``model->width`` x ``model->height``,
with data organized in rows first, columns later.
This means that to access pixel ``(x,y)``
one must access the array at position ``x + y*model->width``.

It was mentioned earlier that the ``example`` profile
fills the image by taking the X and Y coordinates
and filling the pixel
with the value ``|(param1 - param2) * param3 * (x - y)|``.
An implementation of this would then look like this:

.. code-block:: c
 :linenos:
 :emphasize-lines: 8,11,15,18-20

 void profit_make_example(profit_profile *profile, profit_model *model, double *image) {

     double x, y;
     unsigned int i, j;
     double half_xbin = model->xbin/2.;
     double half_ybin = model->ybin/2.;

     profit_example_profile *e = (profit_example_profile *)profile;

     x = 0;
     for (i=0; i < model->width; i++) {
         x += half_xbin;

         y = 0;
         for (j=0; j < model->height; j++) {
            y += half_ybin;

            if ( !model->calcmask || model->calcmask[i + j*model->width] ) {
               double val = fabs( (e->param1 - e->param2) * e->param3 * (x - y) );
               image[i + j*model->width] = val;
            }

            y += half_ybin;
         }
         x += half_xbin;
     }
  }

The code above performs the following steps:

#. First of all, on line 8 we cast the profile argument
   into our profile data type
   so we can reference our specific parameters.
#. On line 11 we loop around the X axis.
   ``i`` is the horizontal pixel index on the image
   and spans from 0 to ``model->width``.
   At the same time we keep track of ``x``,
   which is a floating point number representing
   the horizontal image coordinate
   used to evaluate the profile on that pixel.
   See :doc:`coordinates` for more details
   on the coordinate system used by *libprofit*.
#. Similarly, on line 15 we loop around the Y axis.
#. The model might specify a calculation mask,
   indicating that some pixels should not be calculated,
   which is checked in line 18
#. Being now on a given X and Y coordinate,
   we evaluate our profile on line 19.
#. Finally on line 20 we store the evaluated profile
   on the corresponding pixel of the image.

Creation
^^^^^^^^

Last but not least we look at the creation function.
Its signature looks like this::

 profit_profile *profit_create_example(void);

The creation function creates a new profile structure
and populates it with its default values.
It also associates every new created profile
with the two previously mentioned functions.

For this example the code would look like this::

 profit_profile *profit_create_example() {

    profit_example_profile *e = (profit_example_profile*) malloc(sizeof(profit_example_profile));
    e->profile.init_profile = &profit_init_example;
    e->profile.make_profile = &profit_make_example;

    e->param1 = 1.;
    e->param2 = 2.;
    e->param3 = 3;

    return (profit_profile *)e;
 }

The first line of the function
reserves memory for the new profile.
The second and third lines
bind the initialization and evaluation functions
to your new profile.
Then default values are assigned to each of the parameters,
and finally the new profile,
cast as a ``profit_profile``, is returned.

The signature of the creation function
is also the only one of the three
that needs to be put into the profile's ``.h`` file.
This is necessary for :ref:`wiring_up`.

.. _wiring_up:

Wiring up
---------

To finally wire up your new profile with the rest of *libprofit*
you need to give it a name.
This is done at the ``profit.c`` file.
Open it in an editor
and look for the ``_all_profiles`` array.
This array lists the creation functions of all profiles,
and associate them with an name.

To add the ``example`` profile the following line must thus be added,
just before the end of the array and the *sentinel* item::

 {"example", profit_create_example},

In order to be able to "see" the creation function
the ``example.h`` file must also be included,
which is done earlier on in ``profit.c``::

 #include "example.h"

Finally, and because our compilation system is still very basic,
you need to manually add the new ``.c`` file
to the list of files to be compiled.
This is done by adding it to the ``OBJS`` list
in the ``Makefile``::

 OBJS = sersic.o profit.o [...] example.o


Full example
------------

Below are the full new files
that have been described below.
``example.h`` contains the new data type definition,
plus the signature of the creation function,
while ``example.c`` contains the implementation
of the creation, initialization and evaluation
of ``example`` profiles.

.. literalinclude:: example/example.h
   :caption: example.h
   :language: c
   :linenos:

.. literalinclude:: example/example.c
   :caption: example.c
   :language: c
   :linenos:
