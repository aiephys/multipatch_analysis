# coding: utf8
"""
For generating a DB table describing short term dynamics.

"""
from __future__ import print_function, division

import os, logging
import numpy as np
from collections import OrderedDict
from ...util import timestamp_to_datetime
from ...dynamics import generate_pair_dynamics
from .pipeline_module import MultipatchPipelineModule
from .pulse_response import PulseResponsePipelineModule


class DynamicsPipelineModule(MultipatchPipelineModule):
    """Generates dynamics analysis for each pair
    """
    name = 'dynamics'
    dependencies = [PulseResponsePipelineModule]
    table_group = ['dynamics']
    
    @classmethod
    def create_db_entries(cls, job, session):
        logger = logging.getLogger(__name__)
        db = job['database']
        job_id = job['job_id']
        logger.debug("Processing job %s", job_id)

        # Load experiment from DB
        expt = db.experiment_from_ext_id(job_id, session=session)
        for pair in expt.pairs.values():
            if pair.has_synapse is not True:
                continue
            
            dynamics = generate_pair_dynamics(pair, db, session)
            session.add(dynamics)
            logger.debug("Finished dynamics for pair %s", pair)
        session.commit()
                    
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        return session.query(db.Dynamics).filter(db.Dynamics.pair_id==db.Pair.id).filter(db.Pair.experiment_id==db.Experiment.id).filter(db.Experiment.ext_id.in_(job_ids)).all()
