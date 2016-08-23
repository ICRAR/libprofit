API
===

.. contents::

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

Profile class
-------------

.. doxygenclass:: profit::Profile
   :project: libprofit
   :members:

Model class
-----------

.. doxygenclass:: profit::Model
   :project: libprofit
   :members:
