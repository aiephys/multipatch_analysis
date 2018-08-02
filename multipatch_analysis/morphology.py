# coding: utf8
"""
For generating a DB table describing cell morphology.

"""
from __future__ import print_function, division

import os, sys, multiprocessing

from .database import database as db
from .database import TableGroup
from .pipette_metadata import PipetteMetadata
from . import config


class MorphologyTableGroup(TableGroup):
    schemas = {
        'morphology': [
            """Describes morphological properties of cells.
            """,
            ('cell_id', 'cell.id', 'The ID of the cell described by each record', {'index': True, 'unique': True}),
            ('pyramidal', 'bool', 'Whether the experimenter labeled this cell as pyramidal', {'index': True}),
            # TODO: import more features from LIMS
        ],
    }

    def create_mappings(self):
        TableGroup.create_mappings(self)
        
        Morphology = self['morphology']
        
        db.Cell.morphology = db.relationship(Morphology, back_populates="cell", cascade="delete", single_parent=True, uselist=False)
        Morphology.cell = db.relationship(db.Cell, back_populates="morphology", single_parent=True)


morphology_tables = MorphologyTableGroup()


def init_tables():
    global Morphology
    morphology_tables.create_tables()
    Morphology = morphology_tables['morphology']


# create tables in database and add global variables for ORM classes
init_tables()



@db.default_session
def update_morphology(limit=0, expts=None, parallel=True, workers=6, raise_exceptions=False, session=None):
    """Update morphology table for all experiments
    """
    if expts is None:
        expts_ready = session.query(db.Experiment.acq_timestamp).join(db.Electrode).join(db.Cell).distinct().all()
        expts_done = session.query(db.Experiment.acq_timestamp).join(db.Electrode).join(db.Cell).join(Morphology).distinct().all()

        print("Skipping %d already complete experiments" % (len(expts_done)))
        experiments = [e for e in expts_ready if e not in set(expts_done)]

        if limit > 0:
            np.random.shuffle(experiments)
            experiments = experiments[:limit]

        jobs = [(record.acq_timestamp, index, len(experiments)) for index, record in enumerate(experiments)]
    else:
        jobs = [(expt, i, len(expts)) for i, expt in enumerate(expts)]

    if parallel:
        pool = multiprocessing.Pool(processes=workers)
        pool.map(import_morphology, jobs)
    else:
        for job in jobs:
            import_morphology(job, raise_exceptions=raise_exceptions)


def import_morphology(job_info, raise_exceptions=False):
    session = db.Session()
    
    try:
        expt_id, index, n_jobs = job_info
        print("Importing morphology (expt_id=%f): %d/%d" % (expt_id, index, n_jobs))

        expt = db.experiment_from_timestamp(expt_id, session=session)
        path = os.path.join(config.synphys_data, expt.storage_path)
        pip_meta = PipetteMetadata(path)

        for cell_id,cell in expt.cells.items():
            # How the experimenter described the morphology
            user_morpho = pip_meta.pipettes[cell.ext_id].get('morphology')

            if user_morpho in (None, ''):
                pyramidal = None
            elif user_morpho == 'pyr':
                pyramidal = True
            else:
                print("Unknown morphology string: %s" % user_morpho)
                pyramidal = None

            results = {
                'pyramidal': pyramidal,
            }

            # Write new record to DB
            morphology = Morphology(cell_id=cell.id, **results)
            session.add(morphology)
        session.commit()
    except:
        session.rollback()
        print("Error in experiment: %f" % expt_id)
        if raise_exceptions:
            raise
        else:
            sys.excepthook(*sys.exc_info())
