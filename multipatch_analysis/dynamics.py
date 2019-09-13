import numpy as np
import scipy.stats
from .database import default_db as db


def sorted_pulse_responses(pr_recs):
    """Return all PulseResponse records for *pair* that were generated in a standard induction/recovery
    pulse train, sorted by clamp mode, induction frequency, recovery delay, and pulse number.
    
    Returns
    -------
    sorted_recs : dict
        Nested dictionary structure with keys [(clamp_mode, ind_freq, recovery_delay)][recording][pulse_number]::
            {('ic', 50.0, 250e-3): {
                rec1: {1:pr, 2:pr, 3:pr, ...}, 
                rec2: {1:pr, 2:pr, 3:pr, ...}, 
                ...},
            ...}
         
    Examples
    --------
         
    Get all current clamp pulse responses that are the first pulse in a 50Hz, 250ms delay pulse train:
    
        pr_recs = pulse_response_query(pair).all()
        sorted = sorted_pulse_responses(pr_recs)
        prs = [d[1] for d in sorted['ic', 50, 250e-3].values()]
    """

    # group records by (clamp mode, ind_freq, rec_delay), recording, and then by pulse number
    sorted_recs = {}
    for rec in pr_recs:
        stim_key = (rec.patch_clamp_recording.clamp_mode, rec.multi_patch_probe.induction_frequency, rec.multi_patch_probe.recovery_delay)
        sorted_recs.setdefault(stim_key, {})
        sorted_recs[stim_key].setdefault(rec.recording, {})
        pn = rec.stim_pulse.pulse_number
        sorted_recs[stim_key][rec.recording][pn] = rec
        
    return sorted_recs


def pulse_response_query(pair, qc_pass=False, clamp_mode=None, data=False, session=None):
    if session is None:
        session = db.session()
    q = session.query(db.PulseResponse, db.PulseResponseFit, db.StimPulse, db.Recording, db.PatchClampRecording, db.MultiPatchProbe, db.Synapse)
    q = q.join(db.PulseResponseFit, db.PulseResponse.pulse_response_fit)
    q = q.join(db.StimPulse, db.PulseResponse.stim_pulse)
    q = q.join(db.Recording, db.PulseResponse.recording)
    q = q.join(db.PatchClampRecording, db.PatchClampRecording.recording_id==db.Recording.id)
    q = q.join(db.MultiPatchProbe, db.MultiPatchProbe.patch_clamp_recording_id==db.PatchClampRecording.id)
    q = q.join(db.Synapse, db.Synapse.pair_id==db.PulseResponse.pair_id)
    q = q.filter(db.PulseResponse.pair_id==pair.id)

    if data is True:
        q = q.add_column(db.PulseResponse.data)
    
    if clamp_mode is not None:
        q = q.filter(db.PatchClampRecording.clamp_mode==clamp_mode)
    if qc_pass:
        syn_type = pair.synapse.synapse_type
        pr_qc_field = getattr(db.PulseResponse, syn_type + "_qc_pass")
        q = q.filter(pr_qc_field==True)
    
    return q


def generate_pair_dynamics(pair, db, session):
    """Generate a Dynamics table entry for the given pair.
    """
    syn_type = pair.synapse.synapse_type
    
    # load all IC pulse response amplitudes to determine the maximum that will be used for normalization
    pr_query = pulse_response_query(pair, qc_pass=True, clamp_mode='ic')
    pr_recs = pr_query.all()
    # cull out all PRs that didn't get a fit
    pr_recs = [pr_rec for pr_rec in pr_recs if pr_rec.pulse_response_fit.fit_amp is not None]
    
    percentile = 90 if syn_type == 'ex' else 10
    amps = [rec.pulse_response_fit.fit_amp for rec in pr_recs]
    amp_90p = scipy.stats.scoreatpercentile(amps, percentile)

    # load all baseline amplitudes to determine the noise level
    noise_amps = [rec.pulse_response_fit.baseline_fit_amp for rec in pr_recs if rec.pulse_response_fit.baseline_fit_amp is not None]
    noise_90p = scipy.stats.scoreatpercentile(noise_amps, percentile)

    # start new DB record
    dynamics = db.Dynamics(
        pair_id=pair.id,
        pulse_amp_90th_percentile=amp_90p,
        noise_amp_90th_percentile=noise_90p,
    )

    # sort all PRs by recording and stimulus parameters
    sorted_prs = sorted_pulse_responses(pr_recs)

    # calculate 50Hz paired pulse and induction metrics
    metrics = {'stp_initial_50hz': [], 'stp_induction_50hz': [], 'stp_recovery_250ms': []}
    paired_pulse_ratio = []
    
    for key,recs in sorted_prs.items():
        clamp_mode, ind_freq, rec_delay = key
        if ind_freq != 50:
            continue
        for recording, pulses in recs.items():
            if 1 not in pulses or 2 not in pulses:
                continue
            amps = {k:r.pulse_response_fit.fit_amp for k,r in pulses.items()}
            metrics['stp_initial_50hz'].append((amps[2] - amps[1]) / amp_90p)
            if amps[1] != 0:
                paired_pulse_ratio.append(amps[2] / amps[1])

            if any([k not in pulses for k in [1,6,7,8]]):
                continue
            metrics['stp_induction_50hz'].append((np.mean([amps[6], amps[7], amps[8]]) - amps[1]) / amp_90p)
            
    # PPR is a bit put of place here, but we're including it since it's a popular metric used
    # in the literature.
    dynamics.paired_pulse_ratio_50hz = scipy.stats.gmean(paired_pulse_ratio)
    
    # calculate recovery at 250 ms
    for key,recs in sorted_prs.items():
        clamp_mode, ind_freq, rec_delay = key
        if abs(rec_delay - 250e-3) > 5e-3:
            continue
        for recording, pulses in recs.items():
            if any([k not in pulses for k in range(1,13)]):
                continue
            amps = {k:r.pulse_response_fit.fit_amp for k,r in pulses.items()}
            r = [amps[i+8] - amps[i] for i in range(1,5)]
            metrics['stp_recovery_250ms'].append(np.mean(r) / amp_90p)

    for k,v in metrics.items():
        setattr(dynamics, k, np.mean(v))
        setattr(dynamics, k+'_n', len(v))
        setattr(dynamics, k+'_std', np.std(v))

    return dynamics
