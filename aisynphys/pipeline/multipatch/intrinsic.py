# coding: utf8
"""
For generating a table that describes cell intrinisic properties

"""
from __future__ import print_function, division

import traceback, sys
import numpy as np   

from neuroanalysis.util.optional_import import optional_import
extractors_for_sweeps = optional_import('ipfx.data_set_features', 'extractors_for_sweeps')
LongSquareAnalysis = optional_import('ipfx.stimulus_protocol_analysis', 'LongSquareAnalysis')
Sweep, SweepSet = optional_import('ipfx.sweep', ['Sweep', 'SweepSet'])
extract_chirp_features = optional_import('ipfx.chirp_features', 'extract_chirp_features')
get_complete_long_square_features = optional_import('ipfx.bin.features_from_output_json', 'get_complete_long_square_features')

from .pipeline_module import MultipatchPipelineModule
from .experiment import ExperimentPipelineModule
from .dataset import DatasetPipelineModule
from ...nwb_recordings import get_intrinsic_recording_dict, get_pulse_times


SPIKE_FEATURES = [
    'upstroke_downstroke_ratio',
    'upstroke',
    'downstroke',
    'width',
    'peak_v',
    'threshold_v',
    'fast_trough_v',
]
class IntrinsicPipelineModule(MultipatchPipelineModule):
    
    name = 'intrinsic'
    dependencies = [ExperimentPipelineModule, DatasetPipelineModule]
    table_group = ['intrinsic']

    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        job_id = job['job_id']

        # Load experiment from DB
        expt = db.experiment_from_timestamp(job_id, session=session)
        nwb = expt.data
        if nwb is None:
            raise Exception('No NWB data for this experiment')

        n_cells = len(expt.cell_list)
        errors = []
        for cell in expt.cell_list:
            dev_id = cell.electrode.device_id
            recording_dict = get_intrinsic_recording_dict(expt, dev_id)
            lp_recs = recording_dict['If_Curve'] + recording_dict['TargetV']
            results = {}
            lp_results, error = IntrinsicPipelineModule.get_long_square_features(lp_recs, cell_id=cell.ext_id)
            errors += error
            chirp_results, error = IntrinsicPipelineModule.get_chirp_features(recording_dict['Chirp'], cell_id=cell.ext_id)
            errors += error
            # Write new record to DB
            conn = db.Intrinsic(cell_id=cell.id, **lp_results, **chirp_results)
            session.add(conn)

        return errors

    @staticmethod
    def get_chirp_features(recordings, cell_id=''):
        errors = []
        if len(recordings) == 0:
            errors.append('No chirp sweeps for cell %s' % cell_id)
            return {}, errors
             
        sweep_list = []
        for rec in recordings:
            sweep = MPSweep(rec)
            if sweep is not None:
                sweep_list.append(sweep)
        
        if len(sweep_list) == 0:
            errors.append('No sweeps passed qc for cell %s' % cell_id)
            return {}, errors

        sweep_set = SweepSet(sweep_list) 
        try:
            all_chirp_features = extract_chirp_features(sweep_set)
            results = {
                'chirp_peak_freq': all_chirp_features['peak_freq'],
                'chirp_3db_freq': all_chirp_features['3db_freq'],
                'chirp_peak_ratio': all_chirp_features['peak_ratio'],
            }
        except Exception as exc:
            errors.append('Error processing chirps for cell %s: %s' % (cell_id, str(exc)))
            results = {}
        
        return results, errors

    @staticmethod
    def get_long_square_features(recordings, cell_id=''):
        errors = []
        if len(recordings) == 0:
            errors.append('No long pulse sweeps for cell %s' % cell_id)
            return {}, errors

        min_pulse_dur = np.inf
        sweep_list = []
        for rec in recordings:
            pulse_times = get_pulse_times(rec)
            if pulse_times is None:
                continue
            
            # pulses may have different durations as well, so we just use the smallest duration
            start, end = pulse_times
            min_pulse_dur = min(min_pulse_dur, end-start)
            
            sweep = MPSweep(rec, -start)
            if sweep is not None:
                sweep_list.append(sweep)
        
        if len(sweep_list) == 0:
            errors.append('No sweeps passed qc for cell %s' % cell_id)
            return {}, errors

        sweep_set = SweepSet(sweep_list)    
        spx, spfx = extractors_for_sweeps(sweep_set, start=0, end=min_pulse_dur)
        lsa = LongSquareAnalysis(spx, spfx, subthresh_min_amp=-200)
        
        try:
            analysis = lsa.analyze(sweep_set)
        except Exception as exc:
            errors.append('Error running IPFX analysis for cell %s: %s' % (cell_id, str(exc)))
            return {}, errors
        
        analysis_dict = lsa.as_dict(analysis)
        output = get_complete_long_square_features(analysis_dict) 
        
        results = {
            'rheobase': output['rheobase_i'] * 1e-9, #unscale from pA,
            'fi_slope': output['fi_fit_slope'] * 1e-9, #unscale from pA,
            'input_resistance': output['input_resistance'] * 1e6, #unscale from MOhm,
            'sag': output['sag'],
            'adaptation_index': output['adapt_mean'],
            'upstroke_downstroke_ratio': output['upstroke_downstroke_ratio_hero'],
            'upstroke': output['upstroke_hero'],
            'downstroke': output['downstroke_hero'],
            'width': output['width_hero'],
            'threshold_v': output['threshold_v_hero'],

            'peak_deltav': output['peak_deltav_hero'],
            'fast_trough_deltav': output['fast_trough_deltav_hero'],

            'isi_adapt_ratio': output['isi_adapt_ratio_hero'],
            'upstroke_adapt_ratio': output['upstroke_adapt_ratio_hero'],
            'downstroke_adapt_ratio': output['downstroke_adapt_ratio_hero'],
            'width_adapt_ratio': output['width_adapt_ratio_hero'],
            'threshold_v_adapt_ratio': output['threshold_v_adapt_ratio_hero'],
        }
        return results, errors

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
    def __init__(self, rec, t0=0):
        # pulses may have different start times, so we shift time values to make all pulses start at t=0
        pri = rec['primary'].copy(t0=t0)
        cmd = rec['command'].copy()
        t = pri.time_values
        v = pri.data * 1e3  # convert to mV
        holding = [i for i in rec.stimulus.items if i.description=='holding current']
        if len(holding) == 0:
            # TODO: maybe log this error
            return None
        holding = holding[0].amplitude
        i = (cmd.data - holding) * 1e12   # convert to pA with holding current removed
        srate = pri.sample_rate
        sweep_num = rec.parent.key
        clamp_mode = rec.clamp_mode  # this will be 'ic' or 'vc'; not sure if that's right

        Sweep.__init__(self, t, v, i, clamp_mode, srate, sweep_number=sweep_num)