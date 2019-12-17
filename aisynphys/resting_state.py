# coding: utf8
from __future__ import print_function, division

from .fitting import fit_avg_pulse_response
from .data import PulseResponseList
from .database import default_db as db
from .util import datetime_to_timestamp 


def resting_state_response_fits(pair, rest_duration):
    """Return curve fits to average pulse responses from *pair* for pulses that are at "resting state",
    meaning that each stimulus included in the average is preceded by a certain minimum period
    with no other presynaptic stimuli.
    
    Parameters
    ----------
    pair : Pair instance
        The cell pair from which to query for pulse response data
    rest_duration : float
        Duration (seconds) of the time window that must be quiescent in order to consider
        the synapse at "resting state".
        
    Returns
    -------
    result : dict
        Dictionary containing averages and fit results for current clamp and voltage clamp
    """

    latency = pair.synapse.latency
    if latency is None:
        return None
    latency_window = [latency - 100e-6, latency + 100e-6]
    syn_typ = pair.synapse.synapse_type
        
    fit_signs = {
        ('ex', 'ic'): 1,
        ('ex', 'vc'): -1,
        ('in', 'ic'): -1,
        ('in', 'vc'): 1,
    }
        
    # 1. Select qc-passed "resting state" PRs
    rest_prs = get_resting_state_responses(pair, rest_duration)
    
    # 2. Average and fit
    fits = {}
    for mode, pr_list in rest_prs.items():
        sign = fit_signs[syn_typ, mode]
        if mode == 'ic':
            init_params = {'rise_time': pair.synapse.psp_rise_time, 'decay_tau': pair.synapse.psp_decay_tau}
        else:
            init_params = {'rise_time': pair.synapse.psc_rise_time, 'decay_tau': pair.synapse.psc_decay_tau}
        init_params = {k:v for k,v in init_params.items() if v is not None}

        if len(pr_list) == 0:
            fit, avg = None, None
        else:
            fit, avg = fit_avg_pulse_response(pr_list, latency_window=latency_window, sign=sign, init_params=init_params)

        fits[mode] = {
            'fit': fit,
            'average': avg,
            'responses': pr_list,
        }
    
    return fits


def get_resting_state_responses(pair, rest_duration, response_duration=10e-3):
    """Return {'ic': PulseResponseList(), 'vc': PulseResponseList()} containing
    all qc-passed, resting-state pulse responses for *pair*.
    """
    syn_typ = pair.synapse.synapse_type

    q = db.query(db.PulseResponse, db.StimPulse, db.Recording, db.PatchClampRecording, db.PulseResponse.data)
    q = q.join(db.StimPulse, db.PulseResponse.stim_pulse)
    q = q.join(db.Recording, db.PulseResponse.recording)
    q = q.join(db.PatchClampRecording, db.PatchClampRecording.recording_id==db.Recording.id)
    q = q.filter(db.PulseResponse.pair_id==pair.id)
    q = q.order_by(db.Recording.start_time, db.StimPulse.onset_time)
    recs = q.all()
    
    rest_prs = {'ic': [], 'vc': []}
    stim_times = [datetime_to_timestamp(rec.Recording.start_time) + rec.StimPulse.onset_time for rec in recs]
    for i,rec in enumerate(recs):
        stim_time = stim_times[i]
        if i > 0:
            last_stim_time = stim_times[i-1]
            if stim_time - last_stim_time < rest_duration:
                # not a resting state pulse, skip
                continue
        if i < len(recs) - 1:
            next_stim_time = stim_times[i+1]
            if next_stim_time - stim_time < response_duration:
                # not enough response time; skip
                continue
        
        qc_pass = getattr(rec.PulseResponse, syn_typ + '_qc_pass')
        if qc_pass is not True:
            continue
        
        rest_prs[rec.PatchClampRecording.clamp_mode].append(rec.PulseResponse)
        
        assert rec.PulseResponse.stim_pulse.pulse_number == 1

    return {k:PulseResponseList(v) for k,v in rest_prs.items()}

