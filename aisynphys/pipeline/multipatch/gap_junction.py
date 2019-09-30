# coding: utf8
from __future__ import print_function, division

import os
import pyqtgraph as pg
from ... import config
from ..pipeline_module import DatabasePipelineModule
from .experiment import ExperimentPipelineModule
from .dataset import DatasetPipelineModule


class GapJunctionPipelineModule(DatabasePipelineModule):
    """Analyze gap junction presence and strength for all pairs per experiment
    """
    name = 'gap_junction'
    dependencies = [ExperimentPipelineModule, DatasetPipelineModule]
    table_group = ['gap_junction']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        expt_id = job['job_id']
        
        expt = db.experiment_from_timestamp(expt_id, session=session)

        for pair in expt.pair_list:
            results = {'gap_junction': False, 'strength': 0}

            # Write new record to DB
            conn = db.GapJunction(pair_id=pair.id, **results)
            session.add(conn)
        
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        q = session.query(db.GapJunction)
        q = q.filter(db.GapJunction.pair_id==db.Pair.id)
        q = q.filter(db.Pair.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.ext_id.in_(job_ids))
        return q.all()
