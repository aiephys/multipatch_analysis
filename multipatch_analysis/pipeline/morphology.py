# coding: utf8
"""
For generating a DB table describing cell morphology.

"""
from __future__ import print_function, division

import os
from collections import OrderedDict
from ..util import timestamp_to_datetime
from .. import database as db
from ..pipette_metadata import PipetteMetadata
from .. import config
from .pipeline_module import DatabasePipelineModule
from .experiment import ExperimentPipelineModule


class MorphologyPipelineModule(DatabasePipelineModule):
    """Imports cell morphology data for each experiment
    """
    name = 'morphology'
    dependencies = [ExperimentPipelineModule]
    table_group = db.morphology_tables
    
    @classmethod
    def create_db_entries(cls, job_id, session):
        
        # Load experiment from DB
        expt = db.experiment_from_timestamp(job_id, session=session)
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
            morphology = db.Morphology(cell_id=cell.id, **results)
            session.add(morphology)
        
    @classmethod
    def job_records(cls, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        return session.query(db.Morphology).filter(db.Morphology.cell_id==db.Cell.id).filter(db.Cell.experiment_id==db.Experiment.id).filter(db.Experiment.acq_timestamp.in_(job_ids)).all()

    @classmethod
    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        # All experiments and their creation times in the DB
        expts = ExperimentPipelineModule.finished_jobs()
        
        # Look up nwb file locations for all experiments
        session = db.Session()
        expt_recs = session.query(db.Experiment.acq_timestamp, db.Experiment.storage_path).all()
        expt_paths = {rec.acq_timestamp: rec for rec in expt_recs}
        session.rollback()
        
        # Return the greater of NWB mod time and experiment DB record mtime
        ready = OrderedDict()
        for expt_id, (expt_mtime, success) in expts.items():
            if success is not True:
                continue
            rec = expt_paths[expt_id]
            pip_file = os.path.join(config.synphys_data, rec.storage_path, 'pipettes.yml')
            pip_mtime = timestamp_to_datetime(os.stat(pip_file).st_mtime)
            ready[rec.acq_timestamp] = max(expt_mtime, pip_mtime)
        return ready
