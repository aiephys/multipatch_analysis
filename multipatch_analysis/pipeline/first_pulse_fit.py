# coding: utf8
from __future__ import print_function, division

import os
from .. import database as db
from .. import config
from .pipeline_module import DatabasePipelineModule
from .connection_strength import ConnectionStrengthPipelineModule
from ..fit_average_first_pulse import fit_first_pulses


class FirstPulseFitPipelineModule(DatabasePipelineModule):
    """Analyze synaptic connection strength for all pairs per experiment
    """
    name = 'first_pulse_fit'
    dependencies = [ConnectionStrengthPipelineModule]
    table_group = db.first_pulse_fit_tables
    
    @classmethod
    def create_db_entries(cls, expt_id, session):
    
        expt = db.experiment_from_timestamp(expt_id, session=session)
        
        for (pre_cell_id, post_cell_id), pair in expt.pairs.items():
            result = fit_first_pulses(pair, pre_cell_id, post_cell_id)
            afpf = db.AvgFirstPulseFit(pair=pair, **result)
            session.add(afpf)
        
    @classmethod
    def job_query(cls, job_ids, session):
        """Return a query that returns records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        q = session.query(db.AvgFirstPulseFit)
        q = q.filter(db.AvgFirstPulseFit.pair_id==db.Pair.id)
        q = q.filter(db.Pair.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.acq_timestamp.in_(job_ids))
        return q
