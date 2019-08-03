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
from ...avg_response_fit import get_pair_avg_fits


class AvgResponseFitPipelineModule(DatabasePipelineModule):
    """Generate fit to response average for all pairs per experiment
    """
    name = 'avg_response_fit'
    dependencies = [ExperimentPipelineModule, DatasetPipelineModule]
    table_group = ['avg_response_fit']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        expt_id = job['job_id']
        
        expt = db.experiment_from_timestamp(expt_id, session=session)

        for pair in expt.pair_list:
            fits = get_pair_avg_fits(pair, session)
            for (mode, holding), fit in fits.items():
                if fit is None:
                    continue
                rec = db.AvgResponseFit(
                    pair_id=pair.id,
                    clamp_mode=mode,
                    holding=holding,
                    nrmse=fit['fit_result'].nrmse(),
                    initial_xoffset=fit['initial_latency'],
                    manual_qc_pass=fit['fit_qc_pass'],
                )

                for k in ['xoffset', 'yoffset', 'amp', 'rise_time', 'decay_tau', 'exp_amp', 'exp_tau']:
                    setattr(rec, 'fit_'+k, fit['fit_result'].best_values[k])

                session.add(rec)
        
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        q = session.query(db.AvgResponseFit)
        q = q.filter(db.AvgResponseFit.pair_id==db.Pair.id)
        q = q.filter(db.Pair.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.ext_id.in_(job_ids))
        return q.all()
