# coding: utf8
from __future__ import print_function, division

import os
import numpy as np
import pyqtgraph as pg
from collections import OrderedDict
from ... import config
from .pipeline_module import MultipatchPipelineModule
from .dataset import DatasetPipelineModule
from .morphology import MorphologyPipelineModule
from ...synapse import generate_synapse_record
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
                for cell, cell_classes in synaptic_cell_class.items():
                    if len(set(cell_classes)) == 1:
                        # all synaptic projections agree on sign
                        syn_class = cell_classes[0]
                    else:
                        # mismatched synaptic sign
                        syn_class = None
                        
                    # previously generated nonsynaptic cell class -- based only on transgenic markers and morphology
                    cell_class_ns = cell.cell_class_nonsynaptic
                    
                    if cell_class_ns is None or syn_class == cell_class_ns:
                        # if cell class was not called previously, or if the synaptic class
                        # matches the previous nonsynaptic class
                        cell.cell_class = syn_class
                    elif syn_class is None:
                        cell.cell_class = cell_class_ns
                    cell_meta = cell.meta.copy()
                    cell_meta['synaptic_cell_class'] = syn_class
                    cell.meta = cell_meta
                    cell.cell_class, cell.cell_class_nonsynaptic = cell._infer_cell_classes()

            if pair.has_polysynapse:
                errors = generate_synapse_record(pair, db, session, notes_rec, syn='poly', max_ind_freq=50)
                all_errors.extend(errors)

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
