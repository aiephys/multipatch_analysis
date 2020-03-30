# coding: utf8
from __future__ import print_function, division

import os
import pyqtgraph as pg
import numpy as np
import scipy.stats as stats
from ... import config
from .pipeline_module import MultipatchPipelineModule
from .experiment import ExperimentPipelineModule
from .dataset import DatasetPipelineModule

padding = 30e-3
duration = 150e-3

class GapJunctionPipelineModule(MultipatchPipelineModule):
    """Analyze gap junction presence and strength for all pairs per experiment
    """
    name = 'gap_junction'
    dependencies = [ExperimentPipelineModule, DatasetPipelineModule]
    table_group = ['gap_junction']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        expt_id = job['job_id']
        
        expt = db.experiment_from_timestamp(expt_id, session=session)

        for pair in expt.pair_list:
            try:
                nwb = pair.experiment.data
            except:
                continue
            pre_dev = pair.pre_cell.electrode.device_id
            post_dev = pair.post_cell.electrode.device_id
            try:
                sweeps = nwb.contents
            except:
                continue
            long_pulse = get_lp_sweeps(sweeps, pre_dev, post_dev)
            if long_pulse is None:
                continue

            pre_pulse = []
            post_pulse = []
            pre_noise = []
            post_noise = []
            qc_pass = 0
            for sweep in long_pulse:
                pre_rec = sweep[pre_dev]
                post_rec = sweep[post_dev]
                
                db_pre_rec = get_db_recording(expt, pre_rec)
                if db_pre_rec is None or db_pre_rec.patch_clamp_recording.qc_pass is False:
                    continue
                
                db_post_rec = get_db_recording(expt, post_rec)
                if db_post_rec is None or db_post_rec.patch_clamp_recording.qc_pass is False:
                    continue

                pulse_times = get_pulse_times(pre_rec)
                if pulse_times is None:
                    continue

                qc_pass+=1

                pulse_start = pulse_times[0] + padding
                pulse_stop = pulse_start + duration
                pulse_win = [pulse_start, pulse_stop]
                
                base_stop = pulse_times[0]
                base_start = base_stop - duration
                base_win = [base_start, base_stop]
                
                base2_stop = base_start - padding
                base2_start = base2_stop - duration
                base2_win = [base2_start, base2_stop]
                
                pre_pulse = get_chunk_diff(pre_rec, base_win, pulse_win, pre_pulse)
                post_pulse = get_chunk_diff(post_rec, base_win, pulse_win, post_pulse)
                pre_noise = get_chunk_diff(pre_rec, base2_win, base_win, pre_noise)
                post_noise = get_chunk_diff(post_rec, base2_win, base_win, post_noise)
                
            if qc_pass < 2:
                continue

            r_pulse, p_pulse = stats.pearsonr(pre_pulse, post_pulse)
            r_noise, p_noise = stats.pearsonr(pre_noise, post_noise)
            cc_pulse = coupling_coeff(pre_pulse, post_pulse)
            cc_noise = coupling_coeff(pre_noise, post_noise)
            
            results = {
                'corr_coeff_pulse': r_pulse,
                'p_val_pulse': p_pulse,
                'corr_coeff_noise': r_noise, 
                'p_val_noise': p_noise,
                'coupling_coeff_pulse': cc_pulse,
                'coupling_coeff_noise': cc_noise,
                }

            # Write new record to DB
            conn = db.GapJunction(pair_id=pair.id, **results)
            session.add(conn)
        
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        q = session.query(db.GapJunction)
        q = q.filter(db.GapJunction.pair_id==db.Pair.id)
        q = q.filter(db.Pair.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.ext_id.in_(job_ids))
        return q.all()

def get_lp_sweeps(sweeps, pre_dev, post_dev):
    long_pulse= []
    for sweep in sweeps:
        devices = sweep.devices
        if pre_dev not in devices or post_dev not in devices:
            continue
        rec = sweep[pre_dev]
        stim_name = rec.stimulus.description
        target_v = 'TargetV' in stim_name
        if_curve = 'If_Curve' in stim_name
        if target_v or if_curve:
            long_pulse.append(sweep)
    if len(long_pulse) == 0:
        return None
    return long_pulse

def get_pulse_times(recording):
    # returns the start and stop time of the target_v long pulse
    stims = recording.stimulus.items
    # newer experiments have 5 elements, the 4th one is the long pulse of the target_v stimulus
    if len(stims) == 5:
        return stims[3].start_time, stims[3].start_time + stims[3].duration
    # older experiments only have metadata for the test pulse, use the following to parse the 
    # wave itself to find the stimlus
    elif len(stims) == 2:
        # skip over test pulse
        start = stims[1].start_time + stims[1].duration
        cmd = recording['command'].time_slice(start, None)
        # find pulse by parsing command wave
        dif = np.diff((cmd.data != cmd.data[0]).astype(int))
        inds = np.argwhere(dif != 0)[:, 0]
        if len(inds) == 0:
            return None
        return cmd.time_at(inds[0]), cmd.time_at(inds[1])
    else:
        return None
    
def get_chunk_diff(rec, win1, win2, array):
    chunk1 = rec['primary'].time_slice(win1[0], win1[1])
    chunk2 = rec['primary'].time_slice(win2[0], win2[1])
    delta = chunk2.data.mean() - chunk1.data.mean()
    array.append(delta)
    return array     

def get_db_recording(expt, recording):
    trodes = {e.device_id: e.id for e in expt.electrodes}
    trode_id = trodes[recording.device_id]
    srecs = {srec.ext_id: srec for srec in expt.sync_recs}
    try:
        srec = srecs[recording.parent.key]
    except KeyError:
        return None
    recs = {rec.electrode_id:rec for rec in srec.recordings}
    return recs.get(trode_id, None)

def coupling_coeff(a, b):
    a_array= np.asarray(a)
    b_array = np.asarray(b)
    x, _,_,_ = np.linalg.lstsq(a_array[:, None], b_array[:, None])
    return x[0][0]