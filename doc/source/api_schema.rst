.. _api_schema:

Database Schema
===============

The synaptic physiology database is provided as sqlite files that can be queried using many relational database tools. Although it is possible to read this dataset using standard SQL queries, we also provide an sqlalchemy model that implements a richer interface to the database. The API reference below is derived from the sqlalchemy model classes, but also doubles as a description of the relational database schema.


Slice
-----

.. autoclass:: aisynphys.database.schema.Slice


Experiment
----------

.. autoclass:: aisynphys.database.schema.Experiment

Cell
----

.. autoclass:: aisynphys.database.schema.Cell

Pair
----

.. autoclass:: aisynphys.database.schema.Pair

Morphology
----------

.. autoclass:: aisynphys.database.schema.Morphology

Synapse
-------

.. autoclass:: aisynphys.database.schema.Synapse

Dynamics
--------

.. autoclass:: aisynphys.database.schema.Dynamics

RestingStateFit
---------------

.. autoclass:: aisynphys.database.schema.RestingStateFit

SyncRec
-------

.. autoclass:: aisynphys.database.schema.SyncRec

Recording
---------

.. autoclass:: aisynphys.database.schema.Recording

PatchClampRecording
-------------------

.. autoclass:: aisynphys.database.schema.PatchClampRecording

MultiPatchProbe
---------------

.. autoclass:: aisynphys.database.schema.MultiPatchProbe

PulseResponse
-------------

.. autoclass:: aisynphys.database.schema.PulseResponse

StimPulse
---------

.. autoclass:: aisynphys.database.schema.StimPulse

StimSpike
---------

.. autoclass:: aisynphys.database.schema.StimSpike

Baseline
--------

.. autoclass:: aisynphys.database.schema.Baseline


