# coding: utf8
from __future__ import print_function, division

import os
import numpy as np
import pyqtgraph as pg
from collections import OrderedDict
from ... import config
from ..pipeline_module import DatabasePipelineModule
from .experiment import ExperimentPipelineModule
from .dataset import DatasetPipelineModule
import aisynphys.data.data_notes_db as notes_db
from ...avg_response_fit import get_pair_avg_fits


class SynapsePipelineModule(DatabasePipelineModule):
    """Generate fit to response average for all pairs per experiment
    """
    name = 'synapse'
    dependencies = [ExperimentPipelineModule, DatasetPipelineModule]
    table_group = ['synapse', 'avg_response_fit']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        expt_id = job['job_id']
        
        expt = db.experiment_from_ext_id(expt_id, session=session)

        for pair in expt.pair_list:

            # look up synapse type from notes db
            notes_rec = notes_db.get_pair_notes_record(pair.experiment.ext_id, pair.pre_cell.ext_id, pair.post_cell.ext_id)
            if notes_rec is None:
                continue
            
            # update upstream pair record
            pair.has_synapse = notes_rec.notes['synapse_type'] in ('ex', 'in')
            pair.has_electrical = notes_rec.notes['gap_junction']

            # only proceed if we have a synapse here            
            if not pair.has_synapse:
                continue
            
            # fit PSP shape against averaged PSPs/PCSs at -70 and -55 mV
            fits = get_pair_avg_fits(pair, session)
            
            # collect values with which to decide on the "correct" kinetic values to report
            latency_vals = []
            rise_vals = {'ic': [], 'vc': []}
            decay_vals = {'ic': [], 'vc': []}
            
            for (mode, holding), fit in fits.items():
                if fit is None:
                    continue
                    
                if fit['fit_qc_pass']:
                    # user says this is a good fit; write down the kinetic parameters and number of responses that went into the average
                    latency_vals.append((fit['fit_result'].best_values['xoffset'], len(fit['responses']['qc_pass'])))
                    rise_vals[mode].append((fit['fit_result'].best_values['rise_time'], len(fit['responses']['qc_pass'])))
                    decay_vals[mode].append((fit['fit_result'].best_values['decay_tau'], len(fit['responses']['qc_pass'])))
                    
                # record this fit in the avg_response_fit table
                rec = db.AvgResponseFit(
                    pair_id=pair.id,
                    clamp_mode=mode,
                    holding=holding,
                    nrmse=fit['fit_result'].nrmse(),
                    initial_xoffset=fit['initial_latency'],
                    manual_qc_pass=fit['fit_qc_pass'],
                    avg_data=fit['average'].data,
                    avg_data_start_time=fit['average'].t0,
                    n_averaged_responses=len(fit['responses']),
                    avg_baseline_noise=fit['avg_baseline_noise'],
                )
                reasons = fit['fit_qc_pass_reasons']
                if len(reasons) > 0:
                    rec.meta = {'fit_qc_pass_reasons': reasons}

                for k in ['xoffset', 'yoffset', 'amp', 'rise_time', 'decay_tau', 'exp_amp', 'exp_tau']:
                    setattr(rec, 'fit_'+k, fit['fit_result'].best_values[k])

                session.add(rec)
            
            # create a DB record for this synapse
            syn = db.Synapse(
                pair_id=pair.id,
                synapse_type=notes_rec.notes['synapse_type'],
            )
            print("add synapse:", pair, pair.id)

            # compute weighted average of latency values
            lvals = np.array([lv[0] for lv in latency_vals])
            nvals = np.array([lv[1] for lv in latency_vals])
            if nvals.sum() != 0:
                latency = (lvals * nvals).sum() / nvals.sum()
                dist = np.abs(lvals - latency)
                # only set latency if the averaged values agree
                if np.all(dist < 150e-6):
                    syn.latency = latency
            
            # compute weighted averages of kinetic parameters
            for mode, pfx in [('ic', 'psp_'), ('vc', 'psc_')]:
                for param, fit_vals in [('rise_time', rise_vals[mode]), ('decay_tau', decay_vals[mode])]:
                    vals = np.array([v[0] for v in fit_vals])
                    nvals = np.array([v[1] for v in fit_vals])
                    if nvals.sum() == 0:
                        avg = None
                    else:
                        avg = (vals * nvals).sum() / nvals.sum()
                    setattr(syn, pfx+param, avg)
            
            session.add(syn)

            # session.flush()
        
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        
        q = session.query(db.Synapse)
        q = q.filter(db.Synapse.pair_id==db.Pair.id)
        q = q.filter(db.Pair.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.ext_id.in_(job_ids))
        recs = q.all()
        
        q = session.query(db.AvgResponseFit)
        q = q.filter(db.AvgResponseFit.pair_id==db.Pair.id)
        q = q.filter(db.Pair.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.ext_id.in_(job_ids))
        recs.extend(q.all())

        return recs

    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        dataset_module = self.pipeline.get_module('dataset')
        finished_datasets = dataset_module.finished_jobs()

        # find most recent modification time listed for each experiment
        notes_recs = notes_db.db.query(notes_db.PairNotes.expt_id, notes_db.PairNotes.modification_time)
        mod_times = {}
        for rec in notes_recs:
            mod_times[rec.expt_id] = max(rec.modification_time, mod_times.get(rec.expt_id, rec.modification_time))

        # combine update times from pair_notes and finished_datasets
        ready = OrderedDict()
        for job, (mtime, success) in finished_datasets.items():
            if job in mod_times:
                mtime = max(mtime, mod_times[job])
            ready[job] = mtime
            
        return ready
