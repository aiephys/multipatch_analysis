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

        if len(pr_list) == 0:
            fit, avg = None, None
        else:
            fit, avg = fit_avg_pulse_response(pr_list, latency_window=latency_window, sign=sign)

        fits[mode] = {
            'fit': fit,
            'average': avg,
            'responses': pr_list,
        }
    
    return fits


def get_resting_state_responses(pair, rest_duration):
    """Return {'ic': PulseResponseList(), 'vc': PulseResponseList()} containing
    all qc-passed, resting-state pulse responses for *pair*.
    """
    syn_typ = pair.synapse.synapse_type

    q = db.query(db.PulseResponse, db.StimPulse, db.Recording, db.PatchClampRecording, db.PulseResponse.data)
    q = q.join(db.StimPulse, db.PulseResponse.stim_pulse)
    q = q.join(db.Recording, db.PulseResponse.recording)
    q = q.join(db.PatchClampRecording, db.PatchClampRecording.recording_id==db.Recording.id)
    q = q.filter(db.PulseResponse.pair_id==pair.id)
    q.order_by(db.Recording.start_time, db.StimPulse.onset_time)
    recs = q.all()
    
    rest_prs = {'ic': [], 'vc': []}
    stim_time = None
    for i,rec in enumerate(recs):
        last_stim_time = stim_time
        stim_time = datetime_to_timestamp(rec.recording.start_time) + rec.stim_pulse.onset_time
        if last_stim_time is not None and stim_time - last_stim_time < rest_duration:
            # not a resting state pulse, skip
            continue
        
        qc_pass = getattr(rec.pulse_response, syn_typ + '_qc_pass')
        if qc_pass is not True:
            continue
        
        rest_prs[rec.patch_clamp_recording.clamp_mode].append(rec.pulse_response)
        
        assert rec.pulse_response.stim_pulse.pulse_number == 1

    return {k:PulseResponseList(v) for k,v in rest_prs.items()}

