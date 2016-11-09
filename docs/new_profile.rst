Adding a profile
================

.. default-domain:: cpp
.. highlight:: cpp
.. namespace:: profit

.. contents:: Contents
   :local:

This section explains
the steps required to add a new profile to *libprofit*.

In a nutshell, to add a new profile one must:

* Create a new subclass of :class:`Profile`
* Write the mandatory methods
* Associate the new profile with a standard name

In all steps below,
a completely artificial ``example`` profile is being added,
This new profile takes three parameters:
``param1`` and ``param2`` are double numbers,
while ``param3`` is an unsigned integer.
The profile fills the image by taking the X and Y coordinates
and filling the pixel
with the value ``|(param1 - param2) * param3 * (x - y)|``
and requires that all parameters are positive or 0.

The data types used in this example
are described in detail in :doc:`api`.

New Class
---------

The first step to add a new profile is to
define the C++ class that will hold all its information.
Any kind of information can be added to the class,
but it is **required** that the class extends
the base :class:`Profile` class.
The class should be defined
in an ``.h`` file in the ``profit`` directory
so it can be used by others,
and should be part of the ``profit`` namespace.

So far, it should look like this:

.. code-block:: cpp

 class ExampleProfile : public Profile {
 private:
     double param1;
     double param2;
     unsigned int param3;
 };

Methods
-------

Each profile requires a minimum of three methods
that need to be written:

* The constructor,
* A method to validate the profile's values, and
* A method to evaluate it.

The two latter are imposed by the base class,
and must be called :func:`validate <Profile::validate>` and
:func:`evaluate <Profile::evaluate>`.

In addition, to be able to receive parameters given by the user,
the :func:`parameter <Profile::parameter>` methods
must be overwritten.

Parameters
^^^^^^^^^^

To receive parameters given by the user
the new class must overwrite the necessary
:func:`parameter <Profile::parameter>` methods from the parent class.
There are several flavours of this methods,
depending on the parameter data type,
so only the necessary ones are required.

In our example we only have parameters
of type `double` and `unsigned int`
so we only need to overwrite those two methods.
This method must call its parent method
to check if it already set a parameter with that name,
in which case it should short-cut and return `true`;
it then should check the parameter name
against its own parameters,
and return either `true` or `false`
if the parameter was set or not.

In our example, `double` parameters are set like this:

.. code-block:: cpp

 void ExampleProfile::parameter(const std::string &name, double value) {

     if( Profile::parameter(name, value) ) {
        return true;
     }

     if( name == "param1" )      { param1 = value; }
     else if( name == "param2" ) { param2 = value; }
     else {
        return false;
     }

     return true;

 }

Validation
^^^^^^^^^^

After parameters are all set,
*libprofit* will call the validation function.
The validation function's responsibility,
as its name implies,
is to validate the inputs of the profile,
checking that they obey the required minimum to make the operation successful.

In the case of the ``example`` profile it was mentioned
that all parameters must be positive,
so the code must test for that.
If a violation occurs, a :class:`invalid_parameter` exception is thrown.
This exception will prevent the profile (and in fact the whole model)
from being evaluated.

An example implementation would thus look like this:

.. code-block:: cpp

 void ExampleProfile::validate() {

     if ( this->param1 < 0 ) {
        throw invalid_parameter("param1 is negative");
     }
     if ( this->param1 < 0 ) {
        throw invalid_parameter("param2 is negative");
     }
     if ( this->param3 < 0 ) {
        throw invalid_parameter("param3 is negative");
     }

 }

Note also that the base :class:`Profile` class has a reference
to the model this profile is part of.
Having access to the model means
that one can validate profile-specific values
against model-global values as well.
For example, if a new restriction is added stating
that the ``example`` profile can only be run on images
that are bigger than 20 x 20
then the following code could be added:

.. code-block:: cpp

 if ( this->model->width < 20 || this->model->height < 20 ) {
     throw invalid_parameter("can't apply example profile to images less than 20x20");
 }

Finally, if a profile needs no validation at all
a validation function must still be provided
with an empty body.

Evaluation
^^^^^^^^^^

Next, we look to the :func:`evaluate <Profile::evaluate>` method.
Its ``image`` argument
corresponds to the surface where the pixels must be drawn.
Each profile in the model receives a different image surface,
so you will start with a clean slate.
The image is already initialized with zeros,
so if your profile doesn't cover the entire image
there is no need to cover unfilled areas with zeros.
The image is a vector of ``double`` values
representing a 2D image with dimensions ``model->width`` x ``model->height``,
with data organized in rows first, columns later.
This means that to access pixel ``(x,y)``
one must access the vector at position ``x + y*model->width``.

It was mentioned earlier that the ``example`` profile
fills the image by taking the X and Y coordinates
and filling the pixel
with the value ``|(param1 - param2) * param3 * (x - y)|``.
An implementation of this would then look like this:

.. code-block:: cpp
 :linenos:
 :emphasize-lines: 10,14,17-19

 void ExampleProfile::evaluate(std::vector<double> &image) {

     Model *model = this->model;
     double x, y;
     unsigned int i, j;
     double half_xbin = this->model->scale_x/2.;
     double half_ybin = this->model->scale_y/2.;

     x = 0;
     for (i=0; i < model->width; i++) {
         x += half_xbin;

         y = 0;
         for (j=0; j < model->height; j++) {
            y += half_ybin;

            if ( !model->calcmask || model->calcmask[i + j*model->width] ) {
               double val = fabs( (this->param1 - this->param2) * this->param3 * (x - y) );
               image[i + j*model->width] = val;
            }

            y += half_ybin;
         }
         x += half_xbin;
     }
  }

The code above performs the following steps:

#. On line 10 we loop around the X axis.
   ``i`` is the horizontal pixel index on the image
   and spans from 0 to ``model->width``.
   At the same time we keep track of ``x``,
   which is a floating point number representing
   the horizontal image coordinate
   used to evaluate the profile on that pixel.
   See :doc:`coordinates` for more details
   on the coordinate system used by *libprofit*.
#. Similarly, on line 14 we loop around the Y axis.
#. The model might specify a calculation mask,
   indicating that some pixels should not be calculated,
   which is checked in line 17
#. Being now on a given X and Y coordinate,
   we evaluate our profile on line 18.
#. Finally on line 19 we store the evaluated profile
   on the corresponding pixel of the image.

Constructor
^^^^^^^^^^^

Last but not least we look at the constructor.
Its signature looks like this:

.. code-block:: cpp

 ExampleProfile(const Model &model, const std::string &name);

The constructor arguments must be passed down to the parent class.
The constructor is also in charge of populating the profile
with its default values.
For this example the code would look like this:

.. code-block:: cpp

 ExampleProfile::ExampleProfile(const Model &model, const std::string &name) :
     Profile(model, name),
     param1(1.),
     param2(2.),
     param3(3)
 {
   // no-op
 }


.. _wiring_up:

Wiring up
---------

To finally wire up your new profile with the rest of *libprofit*
you need to give it a name.
This is done at the ``profit.cpp`` file.
Open it in an editor
and look for the ``Model::add_profile`` method.
This method creates different profile instances based on the given name.
Add a new ``else if`` statement to create your new profile
imitating what is done for the other ones.

To add the ``example`` profile the following lines should thus be added
to the first ``if/else if`` block:

.. code-block:: cpp

 else if ( profile_name == "example" ) {
     profile = static_cast<Profile *>(new ExampleProfile());
 }

In order to be able to "see" the constructor
the ``example.h`` file must also be included,
which is done earlier on in ``profit.cpp``:

.. code-block:: cpp

 #include "profit/example.h"

Finally,
you need to manually add the new ``.cpp`` file
to the list of files to be compiled.
This is done by adding it to the ``PROFIT_SRC`` list
in the ``CMakeLists.txt`` file::

 set(PROFIT_SRC
   [...]
   src/example.cpp
   [...]
 )


Full example
------------

Below are the full new files
that have been described below.
``example.h`` contains the new data type definition,
plus the signature of the creation function,
while ``example.cpp`` contains the implementation
of the creation, validation and evaluation
of ``example`` profiles.

.. literalinclude:: example/example.h
   :caption: example.h
   :language: cpp
   :linenos:

.. literalinclude:: example/example.cpp
   :caption: example.cpp
   :language: cpp
   :linenos:
