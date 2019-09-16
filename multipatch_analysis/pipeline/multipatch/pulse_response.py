# coding: utf8
from __future__ import print_function, division

import os, random
import pyqtgraph as pg
from ... import config
from ..pipeline_module import DatabasePipelineModule
from .dataset import DatasetPipelineModule
from .synapse import SynapsePipelineModule
from ...pulse_response_strength import baseline_query, response_query, measure_response, analyze_response_strength


class PulseResponsePipelineModule(DatabasePipelineModule):
    """Analyze postsynaptic responses for all presynaptic evoked spikes
    """
    name = 'pulse_response'
    dependencies = [DatasetPipelineModule, SynapsePipelineModule]
    table_group = ['pulse_response_fit', 'pulse_response_strength', 'baseline_response_strength']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        expt_id = job['job_id']

        rq = response_query(session)
        bq = baseline_query(session)

        # select just data for the selected experiment
        rq = rq.join(db.SyncRec).join(db.Experiment).filter(db.Experiment.ext_id==expt_id)
        bq = bq.join(db.SyncRec).join(db.Experiment).filter(db.Experiment.ext_id==expt_id)
        response_recs = rq.all()
        baseline_recs = bq.all()
        
        # match a baseline to each response
        baselines_by_recording = {}
        for b in baseline_recs:
            baselines_by_recording.setdefault(b.recording_id, []).append(b)
        for b in baselines_by_recording.values():
            random.shuffle(b)
        
        baselines = []
        for r in response_recs:
            b = baselines_by_recording.get(r.recording_id, [])
            if len(b) == 0:
                baselines.append(None)
            else:
                baselines.append(b.pop())
        
        # best estimate of response amplitude using known latency for this synapse
        for rec,baseline_rec in zip(response_recs, baselines):
            if not rec.has_synapse:
                continue
            response_fit, baseline_fit = measure_response(rec, baseline_rec)
            new_rec = db.PulseResponseFit(pulse_response_id=rec.response_id)
            for fit, prefix in [(response_fit, 'fit_'), (baseline_fit, 'baseline_fit_')]:
                if fit is None:
                    continue
                for k in ['amp', 'yoffset', 'rise_time', 'decay_tau', 'exp_amp']:
                    setattr(new_rec, prefix+k, fit.best_values[k])
                setattr(new_rec, prefix+'latency', fit.best_values['xoffset'])
                setattr(new_rec, prefix+'nrmse', fit.nrmse())
            session.add(new_rec)
            
            # keepalive; this loop can take a long time
            session.query(db.Slice).count()
        

        # "unbiased" response analysis used to predict connectivity
        _compute_strength('pulse_response', response_recs, session, db)
        _compute_strength('baseline', baseline_recs, session, db)
        
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        
        q = session.query(db.PulseResponseFit)
        q = q.filter(db.PulseResponseFit.pulse_response_id==db.PulseResponse.id)
        q = q.filter(db.PulseResponse.pair_id==db.Pair.id)
        q = q.filter(db.Pair.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.ext_id.in_(job_ids))
        fits = q.all()
        
        q = session.query(db.PulseResponseStrength)
        q = q.filter(db.PulseResponseStrength.pulse_response_id==db.PulseResponse.id)
        q = q.filter(db.PulseResponse.pair_id==db.Pair.id)
        q = q.filter(db.Pair.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.ext_id.in_(job_ids))
        prs = q.all()
        
        q = session.query(db.BaselineResponseStrength)
        q = q.filter(db.BaselineResponseStrength.baseline_id==db.Baseline.id)
        q = q.filter(db.Baseline.recording_id==db.Recording.id)
        q = q.filter(db.Recording.sync_rec_id==db.SyncRec.id)
        q = q.filter(db.SyncRec.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.ext_id.in_(job_ids))
        brs = q.all()
        
        return fits+prs+brs


def _compute_strength(source, recs, session, db):
    """Compute per-pulse-response strength metrics
    """
    rec_type = db.PulseResponseStrength if source == 'pulse_response' else db.BaselineResponseStrength
    for rec in recs:
        new_rec = {'%s_id'%source: rec.response_id}
        result = analyze_response_strength(rec, source)
        # copy a subset of results over to new record
        for k in ['pos_amp', 'neg_amp', 'pos_dec_amp', 'neg_dec_amp', 'pos_dec_latency', 'neg_dec_latency', 'crosstalk']:
            new_rec[k] = result[k]
        session.add(rec_type(**new_rec))

    # Bulk insert is not safe with parallel processes
    # if source == 'pulse_response':
    #     session.bulk_insert_mappings(PulseResponseStrength, new_recs)
    # else:
    #     session.bulk_insert_mappings(BaselineResponseStrength, new_recs)

    # just to collect error messages here in case we have made a mistake:
    session.flush()
