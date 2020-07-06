# coding: utf8
"""
For generating a table that describes cell intrinisic properties

"""
from __future__ import print_function, division

import traceback, sys
import numpy as np   
from ipfx.data_set_features import extractors_for_sweeps
from ipfx.stimulus_protocol_analysis import LongSquareAnalysis
from ipfx.ephys_data_set import Sweep, SweepSet
from .pipeline_module import MultipatchPipelineModule
from .experiment import ExperimentPipelineModule
from .dataset import DatasetPipelineModule
from ...nwb_recordings import get_lp_sweeps, get_pulse_times, get_db_recording


class IntrinsicPipelineModule(MultipatchPipelineModule):
    
    name = 'intrinsic'
    dependencies = [ExperimentPipelineModule, DatasetPipelineModule]
    table_group = ['intrinsic']

    @classmethod
    def create_db_entries(cls, job, session):
        errors = []
        db = job['database']
        job_id = job['job_id']

        # Load experiment from DB
        expt = db.experiment_from_timestamp(job_id, session=session)
        nwb = expt.data
        if nwb is None:
            raise Exception('No NWB data for this experiment')
        sweeps = nwb.contents

        n_cells = len(expt.cell_list)
        ipfx_fail = 0
        for cell in expt.cell_list:
            dev_id = cell.electrode.device_id
            target_v, if_curve = get_lp_sweeps(sweeps, dev_id)
            lp_sweeps = target_v + if_curve
            if len(lp_sweeps) == 0:
                errors.append('No long pulse sweeps for cell %s' % cell.ext_id)
                continue
            recs = [rec[dev_id] for rec in lp_sweeps]
            min_pulse_dur = np.inf
            sweep_list = []
            for rec in recs:
                if rec.clamp_mode != 'ic':
                    continue

                db_rec = get_db_recording(expt, rec)
                if db_rec is None or db_rec.patch_clamp_recording.qc_pass is False:
                    continue

                pulse_times = get_pulse_times(rec)
                if pulse_times is None:
                    continue
                
                # pulses may have different durations as well, so we just use the smallest duration
                start, end = pulse_times
                min_pulse_dur = min(min_pulse_dur, end-start)
                
                sweep = MPSweep(rec, -start)
                if sweep is None:
                    continue
                sweep_list.append(sweep)
            
            if len(sweep_list) == 0:
                errors.append('No sweeps passed qc for cell %s' % cell.ext_id)
                continue

            sweep_set = SweepSet(sweep_list)    
            spx, spfx = extractors_for_sweeps(sweep_set, start=0, end=min_pulse_dur)
            lsa = LongSquareAnalysis(spx, spfx, subthresh_min_amp=-200)
            
            try:
                analysis = lsa.analyze(sweep_set)
            except Exception as exc:
                errors.append('Error running IPFX analysis for cell %s: %s' % (cell.ext_id, str(exc)))
                ipfx_fail += 1
                continue

            if ipfx_fail == n_cells and n_cells > 1:
                raise Exception('All cells failed IPFX analysis')
                continue
            
            spike_features = lsa.mean_features_first_spike(analysis['spikes_set'])
            up_down = spike_features['upstroke_downstroke_ratio']
            rheo = analysis['rheobase_i']
            fi_slope = analysis['fi_fit_slope']
            input_r = analysis['input_resistance']
            sag = analysis['sag']
            avg_rate = np.mean(analysis['spiking_sweeps'].avg_rate)
            adapt = np.mean(analysis['spiking_sweeps'].adapt)
            
            results = {
                'upstroke_downstroke_ratio': up_down,
                'rheobase': rheo,
                'fi_slope': fi_slope,
                'input_resistance': input_r,
                'sag': sag,
                'avg_firing_rate': avg_rate,
                'adaptation_index': adapt,
            }

            # Write new record to DB
            conn = db.Intrinsic(cell_id=cell.id, **results)
            session.add(conn)

        return errors

    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        q = session.query(db.Intrinsic)
        q = q.filter(db.Intrinsic.cell_id==db.Cell.id)
        q = q.filter(db.Cell.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.ext_id.in_(job_ids))
        return q.all()

class MPSweep(Sweep):
    """Adapter for neuroanalysis.Recording => ipfx.Sweep
    """
    def __init__(self, rec, t0):
        # pulses may have different start times, so we shift time values to make all pulses start at t=0
        pri = rec['primary'].copy(t0=t0)
        cmd = rec['command'].copy()
        t = pri.time_values
        v = pri.data * 1e3  # convert to mV
        holding = [i for i in rec.stimulus.items if i.description=='holding current']
        if len(holding) == 0:
            return None
        holding = holding[0].amplitude
        i = (cmd.data - holding) * 1e12   # convert to pA with holding current removed
        srate = pri.sample_rate
        sweep_num = rec.parent.key
        clamp_mode = rec.clamp_mode  # this will be 'ic' or 'vc'; not sure if that's right

        Sweep.__init__(self, t, v, i, clamp_mode, srate, sweep_number=sweep_num)