.. _interactive_tools:

Interactive Tools
=================

This package includes several interactive tools that can be used to explore the Synaptic Physiology dataset. 
These tools require a python installation including PyQt in addition to standard scientific python packages (see :ref:`installation` for more detailed setup instructions).

:ref:`Matrix Analyzer <matrix_analyzer>`
----------------------------------------
The Matrix Analyzer allows you to browse and interact with the Synaptic Physiology Dataset through a graphical user interface

.. image:: images/matrix_main.*


synaptic_dynamics
-----------------

Simple tool for displaying PSPs and curve fits in response to 12-pulse spike trains (see our `experimental stimuli <https://portal.brain-map.org/explore/connectivity/synaptic-physiology/synaptic-physiology-experiment-methods/experimental-stimuli>`_ for more information). Requires the full-size database. Invoke from an anaconda prompt with::

    conda activate aisynphys
    cd path/to/aisynphys
    python tools/synaptic_dynamics.py --db-version=synphys_r1.0_2019-08-29_full.sqlite

| Note: this command will attempt to download the full 175 GB database (if you don't have it already). Make sure you have enough free space on your drive. The download can be aborted at any time and resumed by running the command again.


pair_analysis
-------------

Tool used for displaying averages of stimulus-responses grouped by clamp mode (voltage- or current-clamp) and holding potential (-70 or -55 mV). This is used py the synaptic physiology project internally for manually annotating each pair as being connected by a chemical / electrical synapse. Requires the full-size database. Invoke from an anaconda prompt with::

    conda activate aisynphys
    cd path/to/aisynphys
    python tools/pair_analysis.py --db-version=synphys_r1.0_2019-08-29_full.sqlite

| Note: this command will attempt to download the full 175 GB database (if you don't have it already). Make sure you have enough free space on your drive. The download can be aborted at any time and resumed by running the command again.
