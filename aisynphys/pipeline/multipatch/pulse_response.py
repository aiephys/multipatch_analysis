# coding: utf8
from __future__ import print_function, division

import os, random
import pyqtgraph as pg
from ... import config
from .pipeline_module import MultipatchPipelineModule
from .dataset import DatasetPipelineModule
from .synapse import SynapsePipelineModule
from ...pulse_response_strength import measure_response, measure_deconvolved_response, analyze_response_strength


class PulseResponsePipelineModule(MultipatchPipelineModule):
    """Analyze postsynaptic responses for all presynaptic evoked spikes
    """
    name = 'pulse_response'
    dependencies = [DatasetPipelineModule, SynapsePipelineModule]
    table_group = ['pulse_response_fit', 'pulse_response_strength']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        expt_id = job['job_id']

        expt = db.experiment_from_ext_id(expt_id, session=session)
        
        n_synapses = len([p for p in expt.pair_list if p.has_synapse])

        # select pulse responses for the selected experiment
        # also request pulse response data to speed up access
        rq = db.query(db.PulseResponse, db.PulseResponse.data)
        rq = rq.join(db.Recording).join(db.SyncRec).join(db.Experiment)
        rq = rq.filter(db.Experiment.ext_id==expt_id)
        prs = [rec.PulseResponse for rec in rq.all()]
        print("%s: got %d pulse responses" % (expt_id, len(prs)))
        
        # best estimate of response amplitude using known latency for this synapse
        fits = 0
        for pr in prs:
            if not pr.pair.has_synapse:
                continue
            
            response_fit, baseline_fit = measure_response(pr)
            response_dec_fit, baseline_dec_fit = measure_deconvolved_response(pr)
            if response_fit is None and response_dec_fit is None:
                # print("no response/dec fits")
                continue
            
            new_rec = db.PulseResponseFit(pulse_response_id=pr.id)
            
            # Psp fits
            for fit, prefix in [(response_fit, 'fit_'), (baseline_fit, 'baseline_fit_')]:
                if fit is None:
                    continue
                for k in ['amp', 'yoffset', 'rise_time', 'decay_tau', 'exp_amp']:
                    if k not in fit.best_values:
                        continue
                    setattr(new_rec, prefix+k, fit.best_values[k])
                setattr(new_rec, prefix+'latency', fit.best_values['xoffset'])
                setattr(new_rec, prefix+'nrmse', fit.nrmse())

            # Deconvolved fits
            for fit, prefix in [(response_dec_fit, 'dec_fit_'), (baseline_dec_fit, 'baseline_dec_fit_')]:
                if fit is None:
                    continue
                for k in ['amp', 'yoffset', 'rise_time', 'decay_tau', 'nrmse']:
                    if k not in fit:
                        continue
                    setattr(new_rec, prefix+k, fit[k])
                setattr(new_rec, prefix+'latency', fit['xoffset'])
                setattr(new_rec, prefix+'reconv_amp', fit['reconvolved_amp'])

            session.add(new_rec)
            fits += 1
            # keepalive; this loop can take a long time
            session.query(db.Slice).count()
        
        print("  %s: added %d fit records for %d synapses" % (expt_id, fits, n_synapses))

        # "unbiased" response analysis used to predict connectivity
        for pr in prs:
            bl = pr.baseline
            if bl is None:
                blid = None
                sources = ['pulse_response']
            else:
                blid = bl.id
                sources = ['pulse_response', 'baseline']
            rec = db.PulseResponseStrength(pulse_response_id=pr.id, baseline_id=blid)
            for source in sources:
                result = analyze_response_strength(pr, source)
                if result is None:
                    continue
                # copy a subset of results over to new record
                for k in ['pos_amp', 'neg_amp', 'pos_dec_amp', 'neg_dec_amp', 'pos_dec_latency', 'neg_dec_latency', 'crosstalk']:
                    k1 = k if source == 'pulse_response' else 'baseline_' + k
                    setattr(rec, k1, result[k])
            session.add(rec)

        # just to collect error messages here in case we have made a mistake:
        session.flush()
        
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
        
        return fits + prs

