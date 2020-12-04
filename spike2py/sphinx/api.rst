API reference
=============

.. module:: spike2py

Signal
------

.. autoclass:: Signal
    :members: remove_offset, filter, calibrate, envelop, normalize, interpolate

Trial
-----

.. autoclass:: Trial
    :members: merge, add_sig

SigInfo
----------

A class to provide information about a signal recorded in Spike2.
Also allows for certain pre-processing steps to be performed,
such as filtering, calibrating and offset removal.

.. autoclass:: SigInfo

TrialInfo
---------

.. autoclass:: TrialInfo

FiltInfo
----------

.. autoclass:: FiltInfo


CalibInfo
---------------

.. autoclass:: CalibInfo


OffsetInfo
----------

.. autoclass:: OffsetInfo

NormInfo
----------

.. autoclass:: NormInfo
