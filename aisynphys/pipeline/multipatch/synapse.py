# coding: utf8
from __future__ import print_function, division

import os, logging
import numpy as np
import pyqtgraph as pg
from collections import OrderedDict
from ... import config
from .pipeline_module import MultipatchPipelineModule
from .dataset import DatasetPipelineModule
from .morphology import MorphologyPipelineModule
from aisynphys.avg_response_fit import get_pair_avg_fits
import aisynphys.data.data_notes_db as notes_db


class SynapsePipelineModule(MultipatchPipelineModule):
    """Basic analysis applied to all cell pairs that have a chemical synapse.

    For all cell pairs, this module first records manual synapse and gap junction calls (from notes database)
    in upstream pair.has_synapse and pair.has_electrical.

    If a chemical synapse is present, then collect any qc-passed pulse responses and sort into 4 categories: (ic -70mV),
    (ic -55mV), (vc -70mV), and (vc -55mV). Pulse responses are averaged within each category and curve-fit; these fit
    parameters are recorded along with a qc-pass/fail flag (see aisynphys.avg_response_fit for more on that topic).
    Fit parameters are stored in the avg_response_fit table.

    A record is also added to the synapse table containing the latency and weighted averages of rise_time / decay_tau.
    Only qc-passed fit data are included in these kinetic parameters; in cases with insufficient qc-passed data,
    these values are left empty. The synapse.psp_amplitude and synapse.psc_amplitude parameters are later filled in 
    by the resting_state module.
    """
    name = 'synapse'
    dependencies = [DatasetPipelineModule, MorphologyPipelineModule]
    table_group = ['synapse', 'avg_response_fit', 'poly_synapse']
    
    @classmethod
    def create_db_entries(cls, job, session):
        all_errors = []
        db = job['database']
        expt_id = job['job_id']
        
        expt = db.experiment_from_ext_id(expt_id, session=session)

        # keep track of whether cells look like they should be inhibitory or excitatory based on synaptic projections
        synaptic_cell_class = {}

        for pair in expt.pair_list:
            try:
                # look up synapse type from notes db
                notes_rec = notes_db.get_pair_notes_record(pair.experiment.ext_id, pair.pre_cell.ext_id, pair.post_cell.ext_id)
                if notes_rec is None:
                    continue
                
                # update upstream pair record
                mono_synapse = notes_rec.notes['synapse_type'] in ('ex', 'in')
                pair.has_synapse = mono_synapse
                poly_synapse = notes_rec.notes.get('polysynaptic_type') in ('ex', 'in', 'mix')
                pair.has_polysynapse = poly_synapse
                pair.has_electrical = notes_rec.notes['gap_junction']

                
                if pair.has_synapse:
                    errors = generate_synapse_record(pair, db, session, notes_rec, syn='mono', max_ind_freq=50)
                    
                    all_errors.extend(errors)

                    pre_cell_class = notes_rec.notes['synapse_type']
                    if pre_cell_class is not None:
                        synaptic_cell_class.setdefault(pair.pre_cell, []).append(pre_cell_class)
                
                    # update cell_class if this is a monosynaptic response:
                    # for cell, cell_classes in synaptic_cell_class.items():
                    cell = pair.pre_cell
                    cell_classes = synaptic_cell_class[cell]
                    if len(set(cell_classes)) == 1:
                        # all synaptic projections agree on sign
                        syn_class = cell_classes[0]
                    else:
                        # mismatched synaptic sign
                        syn_class = 'mixed'
                        
                    # previously generated nonsynaptic cell class -- based only on transgenic markers and morphology
                    # cell_class_ns = cell.cell_class_nonsynaptic
                    
                    # if cell_class_ns is None or syn_class == cell_class_ns:
                    #     # if cell class was not called previously, or if the synaptic class
                    #     # matches the previous nonsynaptic class
                    #     cell.cell_class = syn_class
                    # elif syn_class is None:
                    #     cell.cell_class = cell_class_ns
                    cell_meta = cell.meta.copy()
                    cell_meta['synaptic_cell_class'] = syn_class
                    cell.meta = cell_meta
                    cell.cell_class, _ = cell._infer_cell_classes()
                        # if pair.id==85255 and cell.id==15320:
                        #     sadf    
                if pair.has_polysynapse:
                    errors = generate_synapse_record(pair, db, session, notes_rec, syn='poly', max_ind_freq=50)
                    all_errors.extend(errors)
            except Exception:
                print(f"Error processing pair: {pair}")
                raise

        session.commit()
        return all_errors
        
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

        q = session.query(db.PolySynapse)
        q = q.filter(db.PolySynapse.pair_id==db.Pair.id)
        q = q.filter(db.Pair.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.ext_id.in_(job_ids))
        recs.extend(q.all())
        
        q = session.query(db.AvgResponseFit)
        q = q.filter(db.AvgResponseFit.synapse_id==db.Synapse.id)
        q = q.filter(db.AvgResponseFit.poly_synapse_id==db.PolySynapse.id)
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
            ready[job] = {'dep_time': mtime}
            
        return ready


def generate_synapse_record(pair, db, session, notes_rec, syn='mono', max_ind_freq=50):
    """Generate synapse record and associated avg_response_fits.
    'syn' input denotes whether this is a mono- or poly-synaptic event"""
    logger = logging.getLogger(__name__)
    errors = []

    assert syn in ('mono', 'poly'), "Synapse type must be one of 'mono' or 'poly'"
    # create a DB record for this synapse
    if syn == 'mono':
        syn_entry = db.Synapse(
            pair_id=pair.id,
            synapse_type=notes_rec.notes['synapse_type'],
        )
        logger.info("add synapse: %s %s", pair, pair.id)
    if syn == 'poly':
        syn_entry = db.PolySynapse(
            pair_id=pair.id,
            synapse_type=notes_rec.notes['polysynaptic_type'],
        )
        logger.info("add polysynapse: %s %s", pair, pair.id)
    

    # fit PSP shape against averaged PSPs/PCSs at -70 and -55 mV
    #   - selected from <= 50Hz trains
    #   - must pass ex_qc_pass or in_qc_pass
    #   - must have exactly 1 pre spike with onset time
    fits = get_pair_avg_fits(pair, session, max_ind_freq=max_ind_freq)
    fits_decay_20hz = get_pair_avg_fits(pair, session, max_ind_freq=20)
    # This generates a structure like:
    # {(mode, holding): {
    #     'traces': , 
    #     'average', 
    #     'fit_params',
    #     'initial_latency',
    #     'fit_qc_pass',
    #     'expected_fit_params',
    #     'avg_baseline_noise',
    #     }, 
    # }
    
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
            clamp_mode=mode,
            holding=holding,
            nrmse=fit['fit_result'].nrmse(),
            initial_xoffset=fit['initial_latency'],
            manual_qc_pass=fit['fit_qc_pass'],
            avg_data=fit['average'].data,
            avg_data_start_time=fit['average'].t0,
            n_averaged_responses=len(fit['responses']['qc_pass']),
            avg_baseline_noise=fit['avg_baseline_noise'],
            meta={'expected_fit_params': fit['expected_fit_params'], 'expected_fit_pass': fit['expected_fit_pass']},
        )
        if syn == 'mono':
            rec.synapse = syn_entry
        if syn == 'poly':
            rec.poly_synapse = syn_entry
        
        reasons = fit['fit_qc_pass_reasons']
        if len(reasons) > 0:
            rec.meta = {'fit_qc_pass_reasons': reasons}
            errors.append("Fit errors for %s %s %s: %s" % (pair.experiment.ext_id, pair.pre_cell.ext_id, pair.post_cell.ext_id, '\n'.join(reasons)))

        for k in ['xoffset', 'yoffset', 'amp', 'rise_time', 'decay_tau', 'exp_amp', 'exp_tau']:
            setattr(rec, 'fit_'+k, fit['fit_result'].best_values[k])

        # for decay tau in IC mode we also use trains up to 20Hz, note these are not qc'ed
        if mode == 'ic':
            fit_decay = fits_decay_20hz[(mode, holding)]
            if fit_decay is not None:    
                rec.meta['fit_decay_tau_20hz'] = fit_decay['fit_result'].best_values['decay_tau']

        session.add(rec)
    
    # compute weighted average of latency values
    lvals = np.array([lv[0] for lv in latency_vals])
    nvals = np.array([lv[1] for lv in latency_vals])
    if nvals.sum() != 0:
        latency = (lvals * nvals).sum() / nvals.sum()
        dist = np.abs(lvals - latency)
        # only set latency if the averaged values agree
        if np.all(dist < 200e-6):
            syn_entry.latency = latency
        else:
            errors.append("latency mismatch on %s %s %s" % (pair.experiment.ext_id, pair.pre_cell.ext_id, pair.post_cell.ext_id))
    else:
        errors.append("%s %s: No latency values available for this synapse" % (pair.pre_cell.ext_id, pair.post_cell.ext_id))
    
    # compute weighted averages of kinetic parameters
    for mode, pfx in [('ic', 'psp_'), ('vc', 'psc_')]:
        for param, fit_vals in [('rise_time', rise_vals[mode]), ('decay_tau', decay_vals[mode])]:
            vals = np.array([v[0] for v in fit_vals])
            nvals = np.array([v[1] for v in fit_vals])
            if nvals.sum() == 0:
                errors.append("%s %s: No %s %s values available for this synapse" % (pair.pre_cell.ext_id, pair.post_cell.ext_id, mode, param))
                avg = None
            else:
                avg = (vals * nvals).sum() / nvals.sum()
            setattr(syn_entry, pfx+param, avg)
        
    session.add(syn_entry)
    return errors
