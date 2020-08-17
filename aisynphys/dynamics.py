import logging
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
        stim_key = (rec.PatchClampRecording.clamp_mode, rec.MultiPatchProbe.induction_frequency, rec.MultiPatchProbe.recovery_delay)
        sorted_recs.setdefault(stim_key, {})
        sorted_recs[stim_key].setdefault(rec.Recording, {})
        pn = rec.StimPulse.pulse_number
        sorted_recs[stim_key][rec.Recording][pn] = rec
        
    return sorted_recs


def pulse_response_query(pair, qc_pass=False, clamp_mode=None, data=False, spike_data=False, session=None):
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
    q = q.order_by(db.Recording.start_time, db.StimPulse.onset_time)

    if data is True:
        q = q.add_column(db.PulseResponse.data)
        
    if spike_data is True:
        q = q.add_column(db.StimPulse.data.label('spike_data'))
        q = q.add_column(db.StimPulse.data_start_time.label('spike_data_start_time'))
    
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
    logger = logging.getLogger(__name__)
    logger.info('generate dynamics for %s', pair)
    syn_type = pair.synapse.synapse_type
    
    amp_field = 'dec_fit_reconv_amp'
    baseline_amp_field = 'baseline_' + amp_field
    
    # load all IC pulse response amplitudes to determine the maximum that will be used for normalization
    pr_query = pulse_response_query(pair, qc_pass=False, clamp_mode='ic', session=session)
    pr_recs = pr_query.all()
    # cull out all PRs that didn't get a fit or failed qc
    qc_field = syn_type + '_qc_pass'
    passed_pr_recs = [pr_rec for pr_rec in pr_recs if getattr(pr_rec.PulseResponse, qc_field) and getattr(pr_rec.PulseResponseFit, amp_field) is not None]
    
    percentile = 90 if syn_type == 'ex' else 10
    # dec_fit_reconv_amp generally has much lower noise than fit_amp:
    amps = [getattr(rec.PulseResponseFit, amp_field) for rec in passed_pr_recs]
    amp_90p = scipy.stats.scoreatpercentile(amps, percentile)

    # load all baseline amplitudes to determine the noise level
    noise_amps = np.array([getattr(rec.PulseResponseFit, baseline_amp_field) for rec in passed_pr_recs if getattr(rec.PulseResponseFit, baseline_amp_field) is not None])
    noise_std = noise_amps.std()
    noise_90p = scipy.stats.scoreatpercentile(noise_amps, percentile)

    # start new DB record
    dynamics = db.Dynamics(
        pair_id=pair.id,
        pulse_amp_90th_percentile=amp_90p,
        noise_amp_90th_percentile=noise_90p,
        noise_std=noise_std,
    )

    # sort all PRs by recording and stimulus parameters
    #   [(clamp_mode, ind_freq, recovery_delay)][recording][pulse_number]
    sorted_prs = sorted_pulse_responses(passed_pr_recs)

    # calculate 50Hz paired pulse and induction metrics
    metrics = {'stp_initial_50hz': [], 'stp_induction_50hz': [], 'stp_recovery_250ms': []}
    paired_pulse_ratio = []
    all_metrics = {}
    induction = {}
    recovery = {}
    delays = [125e-3, 250e-3, 500e-3, 1000e-3, 2000e-3, 4000e-3]
    
    for key,recs in sorted_prs.items():
        clamp_mode, ind_freq, rec_delay = key
        # if ind_freq != 50:
        #     continue
        if ind_freq not in induction.keys():
            induction[ind_freq] = {'stp_induction': [], 'stp_initial': []}
        for recording, pulses in recs.items():
            if 1 not in pulses or 2 not in pulses:
                continue
            amps = {k:getattr(r.PulseResponseFit, amp_field) for k,r in pulses.items()}
            initial = (amps[2] - amps[1]) / amp_90p
            # we separate out 50Hz into its own column because the most data is here
            if ind_freq == 50:
                metrics['stp_initial_50hz'].append(initial)
                if amps[1] != 0:
                    paired_pulse_ratio.append(amps[2] / amps[1])
            induction[ind_freq]['stp_initial'].append(initial)
            if any([k not in pulses for k in [1,6,7,8]]):
                continue
            ind = (np.mean([amps[6], amps[7], amps[8]]) - amps[1]) / amp_90p
            if ind_freq == 50:
                metrics['stp_induction_50hz'].append(ind)
            induction[ind_freq]['stp_induction'].append(ind)
            
    # PPR is a bit out of place here, but we're including it since it's a popular metric used
    # in the literature.
    dynamics.paired_pulse_ratio_50hz = scipy.stats.gmean(paired_pulse_ratio)
    
    # calculate recovery at 250 ms
    for key,recs in sorted_prs.items():
        clamp_mode, ind_freq, rec_delay = key
        if rec_delay is None:
            continue 
        check_delays = [abs(rec_delay - d) < 5e-3 for d in delays]
        if sum(check_delays) != 1:
            continue
        delay = delays[check_delays[0]]
        if delay not in 
        for recording, pulses in recs.items():
            if any([k not in pulses for k in range(1,13)]):
                continue
            amps = {k:getattr(r.PulseResponseFit, amp_field) for k,r in pulses.items()}
            r = [amps[i+8] - amps[i] for i in range(1,5)]
            if ind_freq == 50 and delay == 250e-3
                metrics['stp_recovery_250ms'].append(np.mean(r) / amp_90p)


    for k,v in metrics.items():
        setattr(dynamics, k, np.mean(v))
        setattr(dynamics, k+'_n', len(v))
        setattr(dynamics, k+'_std', np.std(v))
        
    # Measure PSP variability -- we want a metric something like the coefficient of variation, but 
    # normalized against the 90th% amplitude, and with measurement noise subtracted out. This
    # ends up looking like:
    #     sqrt(amp_stdev^2 - noise_stdev^2) / abs(amp_90th_percentile)
        
    # Variability at resting state:
    resting_amps = []
    for pr_rec in pr_recs:
        if pr_rec.StimPulse.previous_pulse_dt > 8.0 and getattr(pr_rec.PulseResponse, qc_field):
            resting_amps.append(getattr(pr_rec.PulseResponseFit, amp_field))
    
    if len(resting_amps) == 0:
        logger.info("%s: no resting amps; bail out", pair)
        return dynamics
        
    def variability(x):
        return np.log((np.std(x)**2 - noise_std**2)**0.5 / abs(np.mean(x)))

    dynamics.variability_resting_state = variability(resting_amps)

    # Variability in STP-induced state (5th-8th pulses)
    pulse_amps = {
        (2,3): [],
        (5,9): [],
    }
    # collect pulse amplitudes in each category
    for key,recs in sorted_prs.items():
        clamp_mode, ind_freq, rec_delay = key
        if ind_freq != 50:
            continue
        for recording, pulses in recs.items():
            for pulse_n, amps in pulse_amps.items():
                if any([k not in pulses for k in range(1,pulse_n[-1])]):
                    continue
                for n in range(*pulse_n):
                    amps.append(getattr(pulses[n].PulseResponseFit, amp_field))
    
    # normalize
    pulse_var = {n:(variability(a) if len(a) > 0 else np.nan) for n,a in pulse_amps.items()}
        
    # record changes in vairabilty
    dynamics.variability_second_pulse_50hz = pulse_var[2,3]
    dynamics.variability_stp_induced_state_50hz = pulse_var[5,9]
    dynamics.variability_change_initial_50hz = pulse_var[2,3] - dynamics.variability_resting_state
    dynamics.variability_change_induction_50hz = pulse_var[5,9] - dynamics.variability_resting_state
    
    # Look for evidence of vesicle depletion -- correlations between adjacent events in 50Hz pulses 5-8.
    ev1_amp = []
    ev2_amp = []
    pulse_amps = [[] for i in range(9)]
    for key,recs in sorted_prs.items():
        clamp_mode, ind_freq, rec_delay = key
        if ind_freq != 50:
            continue
        for recording, pulses in recs.items():
            for i in range(1,9):
                if i not in pulses:
                    break
                if i < 6:
                    continue
                ev1_amp.append(getattr(pulses[i-1].PulseResponseFit, amp_field))
                ev2_amp.append(getattr(pulses[i].PulseResponseFit, amp_field))

    ev1_amp = np.array(ev1_amp)
    ev2_amp = np.array(ev2_amp)
    ev1_amp -= np.median(ev1_amp)
    ev2_amp -= np.median(ev2_amp)
    
    r,p = scipy.stats.pearsonr(ev1_amp, ev2_amp)
    dynamics.paired_event_correlation_r = r
    dynamics.paired_event_correlation_p = p
    
    return dynamics
