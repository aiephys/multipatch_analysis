# coding: utf8
from __future__ import print_function, division

import os
import pyqtgraph as pg
from collections import OrderedDict
from ... import config
from ..pipeline_module import DatabasePipelineModule
from .experiment import ExperimentPipelineModule
from .dataset import DatasetPipelineModule
from .pulse_response import PulseResponsePipelineModule
from ...avg_response_fit import response_query, sort_responses, fit_avg_response, pair_notes_query


class AvgResponseFitPipelineModule(DatabasePipelineModule):
    """Generate fit to response average for all pairs per experiment
    """
    name = 'avg_response_fit'
    dependencies = [ExperimentPipelineModule, DatasetPipelineModule, PulseResponsePipelineModule]
    table_group = ['avg_response_fit']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        expt_id = job['job_id']
        s2 = notes_db.db.session()
        
        expt = db.experiment_from_timestamp(expt_id, session=session)

        for pair in expt.pair_list:
            # Generate fit results for this pair
            q = response_query(session=session, pair=pair)
            pulse_responses = q.all()
            traces, _ = sort_responses(pulse_responses)
            modes, holdings = traces.items()
            holdings = holdings.keys()
            fit_parameters = OrderedDict()
            q2 = pair_notes_query(session=s2, pair=pair)
            notes = q2.all()
            for clamp_mode in modes:
                for holding in holdings:
                fit_parameters[clamp_mode][holding] = {}  
                    if len(notes) == 0:
                        latency = None
                        sign = 'any'
                    elif len(notes) == 1:
                        latency = notes.notes['fit_parameters']['initial'][clamp_mode][holding]['xoffset']
                        sign = notes.notes['synapse_type'] 
                    else:
                        raise Exception('More than one record for this pair %s %s->%s was found in the Pair Notes database' % (expt_id, pre_cell_id, post_cell_id))

                    fit_parameters[clamp_mode][holding], xoffset, _ = fit_avg_response(traces, clamp_mode, holding, latency, sign)
            
            results = {
            'latency': xoffset,
            'fit_parameters': fit_parameters
            }

            # Write new record to DB
            conn = db.AvgResponseFit(pair_id=pair.id, **results)
            session.add(conn)
        
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        q = session.query(db.AvgResponseFit)
        q = q.filter(db.AvgResponseFit.pair_id==db.Pair.id)
        q = q.filter(db.Pair.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.acq_timestamp.in_(job_ids))
        return q.all()
