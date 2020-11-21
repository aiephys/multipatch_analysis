.. _installation:

Installation
============

Get analysis rolling on your own computer by following these instructions to create a dedicated Anaconda environment. During this process you'll install Anaconda and Git as well as clone repositories for data processing and plotting.

#. Install Anaconda and Git:

    * Download miniconda 3 from `https://repo.continuum.io/miniconda/ <https://repo.continuum.io/miniconda/>`_

        * Pick ``Miniconda3-latest-YourOS-x86-64``

    * Download git from `https://git-scm.com/downloads <https://git-scm.com/downloads>`_

        * Be sure to enable ``git from the command line`` during install
        * We also recommend ``check out as-is, commit as-is`` for Windows users

#. Set up conda environment with aisynphys code:

    * From Anaconda prompt, download aisynphys repository::

        > git clone https://github.com/alleninstitute/aisynphys
        > cd aisynphys

    * Create environment::

        > conda env create --name aisynphys --file desktop-environment.yml
        > conda activate aisynphys

    * Install aisynphys::

        > python setup.py develop
        > cd ..

#. Check out the :ref:`interactive tools <interactive_tools>` you can use!