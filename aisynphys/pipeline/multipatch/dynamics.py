# coding: utf8
"""
For generating a DB table describing short term dynamics.

"""
from __future__ import print_function, division

import os
import numpy as np
from collections import OrderedDict
from ...util import timestamp_to_datetime
from ...dynamics import generate_pair_dynamics
from ..pipeline_module import DatabasePipelineModule
from .pulse_response import PulseResponsePipelineModule
from .synapse_prediction import SynapsePredictionPipelineModule



class DynamicsPipelineModule(DatabasePipelineModule):
    """Generates dynamics analysis for each pair
    """
    name = 'dynamics'
    dependencies = [PulseResponsePipelineModule]
    table_group = ['dynamics']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        job_id = job['job_id']

        delays = [125, 250, 500, 1000, 2000, 4000]
        # Load experiment from DB
        expt = db.experiment_from_timestamp(job_id, session=session)
        for pair in expt.pairs.values():
            if pair.has_synapse is not True:
                continue
            
            dynamics = generate_pair_dynamics(pair, db, session)
            session.add(dynamics)
        
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        return session.query(db.Dynamics).filter(db.Dynamics.pair_id==db.Pair.id).filter(db.Pair.experiment_id==db.Experiment.id).filter(db.Experiment.ext_id.in_(job_ids)).all()
