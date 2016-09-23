API
===

.. contents:: Contents
   :local:

*libprofit* has two main data types:
the model class (:class:`Model`)
and the base profile class (:class:`Profile`).
We introduce the base profile first, then the model.

As a small reference, see the following class diagram.

.. image:: images/types.png

Exceptions
----------

.. doxygenclass:: profit::invalid_parameter
   :project: libprofit

Model class
-----------

.. doxygenclass:: profit::Model
   :members:

Profile classes
---------------

.. doxygenclass:: profit::Profile
   :members:

.. doxygenclass:: profit::RadialProfile
   :members:

.. doxygenclass:: profit::SersicProfile
   :members:

.. doxygenclass:: profit::MoffatProfile
   :members:

.. doxygenclass:: profit::FerrerProfile
   :members:

.. doxygenclass:: profit::CoreSersicProfile
   :members:

.. doxygenclass:: profit::BrokenExponentialProfile
   :members:

.. doxygenclass:: profit::KingProfile
   :members:

.. doxygenclass:: profit::PsfProfile
   :members:

.. doxygenclass:: profit::SkyProfile
   :members:
