# coding: utf8
from __future__ import print_function, division

import os
import pyqtgraph as pg
from ... import config
from ..pipeline_module import DatabasePipelineModule
from .experiment import ExperimentPipelineModule
from .dataset import DatasetPipelineModule
from .pulse_response import PulseResponsePipelineModule
from ...synapse_prediction import get_amps, get_baseline_amps, analyze_pair_connectivity


class SynapsePredictionPipelineModule(DatabasePipelineModule):
    """Analyze synaptic connection strength for all pairs per experiment
    """
    name = 'synapse_prediction'
    dependencies = [ExperimentPipelineModule, DatasetPipelineModule, PulseResponsePipelineModule]
    table_group = ['synapse_prediction']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        expt_id = job['job_id']
        
        expt = db.experiment_from_timestamp(expt_id, session=session)

        for pair in expt.pair_list:
            # Query all pulse amplitude records for each clamp mode
            amps = {}
            for clamp_mode in ('ic', 'vc'):
                clamp_mode_fg = get_amps(session, pair, clamp_mode=clamp_mode, get_data=True)
                clamp_mode_bg = get_baseline_amps(session, pair, amps=clamp_mode_fg, clamp_mode=clamp_mode, get_data=False)
                amps[clamp_mode, 'fg'] = clamp_mode_fg
                amps[clamp_mode, 'bg'] = clamp_mode_bg
            
            # Generate summary results for this pair
            results = analyze_pair_connectivity(amps)
            
            if results is None:
                # no data to analyze
                continue

            # Write new record to DB
            conn = db.SynapsePrediction(pair_id=pair.id, **results)
            session.add(conn)
        
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        q = session.query(db.SynapsePrediction)
        q = q.filter(db.SynapsePrediction.pair_id==db.Pair.id)
        q = q.filter(db.Pair.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.ext_id.in_(job_ids))
        return q.all()
