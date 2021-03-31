# coding: utf8
from __future__ import print_function, division

from .avg_response_fit import fit_avg_pulse_response
from .data import PulseResponseList
from .database import default_db as db
from .util import datetime_to_timestamp 


def resting_state_response_fits(synapse, rest_duration):
    """Return curve fits to average pulse responses from *synapse* for pulses that are at "resting state",
    meaning that each stimulus included in the average is preceded by a certain minimum period
    with no other presynaptic stimuli. 
    
    Parameters
    ----------
    synapse : Synapse or PolySynapse instance
        The synapse from which to query for pulse response data
    rest_duration : float
        Duration (seconds) of the time window that must be quiescent in order to consider
        the synapse at "resting state".
        
    Returns
    -------
    result : dict
        Dictionary containing averages and fit results for current clamp and voltage clamp
    """

    latency = synapse.latency
    
    if latency is None:
        return None
    latency_window = [latency - 100e-6, latency + 100e-6]
    syn_typ = synapse.synapse_type
        
    fit_signs = {
        ('ex', 'ic'): 1,
        ('ex', 'vc'): -1,
        ('in', 'ic'): -1,
        ('in', 'vc'): 1,
    }
        
    # 1. Select qc-passed "resting state" PRs
    rest_prs = get_resting_state_responses(synapse, rest_duration, response_duration=10e-3)
    
    # 2. Average and fit
    fits = {}
    for mode, pr_list in rest_prs.items():
        sign = fit_signs[syn_typ, mode]
        if mode == 'ic':
            init_params = {'rise_time': synapse.psp_rise_time, 'decay_tau': synapse.psp_decay_tau}
        else:
            init_params = {'rise_time': synapse.psc_rise_time, 'decay_tau': synapse.psc_decay_tau}
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


def get_resting_state_responses(synapse, rest_duration, response_duration):
    """Return {'ic': PulseResponseList(), 'vc': PulseResponseList()} containing
    all qc-passed, resting-state pulse responses for *pair*.
    
    The *rest_duration* parameter is used to define the stimuli that count as "resting state":
    any pulse-response that is preceded by a window *rest_duration* seconds long in which there
    are no presynaptic spikes. Typical values here might be a few seconds to tens of seconds to 
    allow the synapse to recover to its resting state.
    """
    syn_typ = synapse.synapse_type
    qc_field = getattr(db.PulseResponse, syn_typ + '_qc_pass')
    
    q = db.query(db.PulseResponse, db.StimPulse, db.Recording, db.PatchClampRecording, db.PulseResponse.data)
    q = q.join(db.StimPulse, db.PulseResponse.stim_pulse)
    q = q.join(db.Recording, db.PulseResponse.recording)
    q = q.join(db.PatchClampRecording, db.PatchClampRecording.recording_id==db.Recording.id)
    q = q.filter(db.PulseResponse.pair_id==synapse.pair_id)
    q = q.filter(db.StimPulse.previous_pulse_dt > rest_duration)
    q = q.filter(qc_field==True)
    q = q.order_by(db.Recording.start_time, db.StimPulse.onset_time)
    recs = q.all()
        
    rest_prs = {'ic': [], 'vc': []}
    for rec in recs:
        pr = rec.PulseResponse
        if pr.post_tseries.duration < response_duration:
            # not enough data; skip
            continue
        rest_prs[rec.PatchClampRecording.clamp_mode].append(pr)

    return {k:PulseResponseList(v) for k,v in rest_prs.items()}
