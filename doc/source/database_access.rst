.. _database_access:

Database Access
===============

The synaptic physiology dataset is implemented as a relational database derived from HDF5 (NWBv1) formatted patch clamp recording data. Although it is possible to access the database using standard SQL queries, we also provide an sqlalchemy model that implements a richer interface to the database. The documentation below provides examples of database access that assume some basic understanding of relational databases and the sqlalchemy ORM.

.. seealso::
    
    * :ref:`Synaptic Physiology Dataset Structure <dataset_structure>`
    * `SQLAlchemy query API <https://docs.sqlalchemy.org/en/13/orm/query.html>`_
    * `SQLAlchemy object relational tutorial <https://docs.sqlalchemy.org/en/13/orm/tutorial.html>`_


Loading a database version
--------------------------

The synaptic physiology database can be loaded using :func:`SynphysDatabase.load_version() <aisynphys.database.SynphysDatabase.load_version>`::

    from aisynphys.database import SynphysDatabase
    db = SynphysDatabase.load_version('synphys_r1.0_2019-08-29_small.sqlite')

This method will take care of downloading and caching the requested database file. If the download is interrupted, it can be automatically resumed by running the function again. The object returned is a :class:`SynphysDatabase <aisynphys.database.SynphysDatabase>` instance that serves as the starting point for all subsequent data access. 


.. note::
    
    * Get a list of all available DB versions with :func:`SynphysDatabase.list_versions() <aisynphys.database.SynphysDatabase.load_version>`
    * Define the location where DB files are downloaded and cached using::

          aisynphys.config.cache_path = "path/to/cache_files"


Querying pairs from the database
--------------------------------

The most common database queries begin by searching for a set of :class:`Pair <aisynphys.database.schema.Pair>` records. The :func:`SynphysDatabase.pair_query <aisynphys.database.SynphysDatabase.pair_query>` method provides a simple interface for constructing this type of query::

    # a query that retrieves all pairs from mouse projects
    query = db.pair_query(project_name=db.mouse_projects)

    # a query that retrieves all pairs from human projects 
    # that are connected by a chemical synapse
    query = db.pair_query(project_name=db.human_projects, synapse=True)

    # a query that retrieves all chemical synapses from mouse projects
    # where the postsynaptic cell is inhibitory
    from aisynphys.cell_class import CellClass
    inhibitory_class = CellClass(cre_type=('pvalb', 'sst', 'vip'))
    query = db.pair_query(project_name=db.mouse_projects, synapse=True, 
                          post_cell_class=inhibitory_class)

The objects returned are `SQLALchemy query objects <https://docs.sqlalchemy.org/en/13/orm/query.html#the-query-object>`_ that can be used to retrieve query results::

    # return a list of all results from the query
    pairs = query.all()

[more on this topic coming soon..]


.. seealso::
    
    * `Example connectivity notebook <https://nbviewer.jupyter.org/github/AllenInstitute/aisynphys/blob/documentation/connectivity.ipynb>`_
    * `Example kinetics notebook <https://nbviewer.jupyter.org/github/AllenInstitute/aisynphys/blob/documentation/synaptic_kinetics.ipynb>`_
    * `Example short-term plasticity notebook <https://nbviewer.jupyter.org/github/AllenInstitute/aisynphys/blob/documentation/short_term_plasticity.ipynb>`_
    

General database queries
------------------------

To construct more general queries, we provide access to SQLAlchemy Session and Mapping classes. Note that all Mapping classes described in the :ref:`database schema <api_schema>` can be accessed as attributes of the :class:`SynphysDatabase <aisynphys.database.SynphysDatabase>` instance::

    # load up a database
    from aisynphys.database import SynphysDatabase
    db = SynphysDatabase.load_version('synphys_r1.0_2019-08-29_small.sqlite')

    # create an SQLAlchemy session
    session = db.session()

    # build a query to select all mouse experiments
    query = session.query(db.Experiment).join(db.Slice).filter(db.Slice.species=='mouse')

    # retrieve the list of experiments
    expts = query.all()

    # from the first experiment returned, print a list of cells and their transgenic type:
    print("Frst experiment returned:", expts[0])
    for cell in expts[0].cells:
        print("Cell %s cre type: %s" % (cell.ext_id, cell.cre_type))

[more on this topic coming soon..]



.. seealso::

    * `Querying in the SQLAlchemy ORM tutorial <https://docs.sqlalchemy.org/en/13/orm/tutorial.html#querying>`_
