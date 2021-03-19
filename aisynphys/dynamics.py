import logging
import numpy as np
import scipy.stats
import pandas as pd
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

    # dec_fit_reconv_amp generally has much lower noise than fit_amp:
    amp_field = 'dec_fit_reconv_amp'
    baseline_amp_field = 'baseline_' + amp_field
    
    # QC criteria
    min_amps_length = 10

    # start new DB record
    dynamics = db.Dynamics(
        pair_id=pair.id,
        meta={},
    )

    syn_type = pair.synapse.synapse_type
    
    # load all IC pulse response amplitudes to determine the maximum that will be used for normalization
    pr_query = pulse_response_query(pair, qc_pass=False, clamp_mode='ic', session=session)
    pr_recs = pr_query.all()
    # cull out all PRs that didn't get a fit or failed qc
    qc_field = syn_type + '_qc_pass'
    passed_pr_recs = [pr_rec for pr_rec in pr_recs if getattr(pr_rec.PulseResponse, qc_field) and getattr(pr_rec.PulseResponseFit, amp_field) is not None]
    
    percentile = 90 if syn_type == 'ex' else 10
    amps = [getattr(rec.PulseResponseFit, amp_field) for rec in passed_pr_recs]
    amp_90p = scipy.stats.scoreatpercentile(amps, percentile)

    # fail QC if there are not enough events (but continue anyway)
    qc_pass = len(amps) >= min_amps_length
    dynamics.qc_pass = qc_pass
    dynamics.n_source_events = len(amps)

    # load all baseline amplitudes to determine the noise level
    noise_amps = np.array([getattr(rec.PulseResponseFit, baseline_amp_field) for rec in passed_pr_recs if getattr(rec.PulseResponseFit, baseline_amp_field) is not None])
    noise_std = noise_amps.std()
    noise_90p = scipy.stats.scoreatpercentile(noise_amps, percentile)

    # fill in new fields
    dynamics.pulse_amp_90th_percentile = amp_90p
    dynamics.noise_amp_90th_percentile = noise_90p
    dynamics.noise_std = noise_std

    # sort all PRs by recording and stimulus parameters
    #   [(clamp_mode, ind_freq, recovery_delay)][recording][pulse_number]
    sorted_prs = sorted_pulse_responses(passed_pr_recs)

    # calculate 50Hz paired pulse and induction metrics for their own column
    col_metrics = {'stp_initial_50hz': [], 'stp_induction_50hz': [], 'stp_recovery_250ms': [], 'stp_recovery_single_250ms': []}
    paired_pulse_ratio = []
    # caclulate dynamics for all frequencie and recovery delays
    #   [(clamp_mode, ind_freq, recovery_delay), {'stp_induction':(mean, std, n),
    #       'stp_initial': (mean, std, n), 'stp_recovery': (mean, std, n)}), ...]
    all_metrics = []
    delays = [125e-3, 250e-3, 500e-3, 1000e-3, 2000e-3, 4000e-3]
    
    for key,recs in sorted_prs.items():
        clamp_mode, ind_freq, rec_delay = key
        # check for a rec_delay and match it to closest known interval in delays
        if rec_delay is None:
            delay = None
        else:
            check_delays = np.array([abs(rec_delay-d) <= 5e-3 for d in delays]) # allow measured delay to be off by 5ms
            if sum(check_delays) == 1:
                delay = np.array(delays)[check_delays][0]
            else:
                delay = rec_delay
        meta = (clamp_mode, ind_freq, delay)
        
        collect_initial = []
        collect_induction = []
        collect_recovery = []
        collect_recovery_single = []
        collect_pulse_amps = [[]] * 12

        for recording, pulses in recs.items():
            # get all pulse amps
            amps = {k:getattr(r.PulseResponseFit, amp_field) for k,r in pulses.items()}

            # calculate metrics if the proper conditions are met
            if 1 in pulses and 2 in pulses:
                initial = (amps[2] - amps[1]) / amp_90p
                collect_initial.append(initial)
                # we separate out 50Hz into its own column because the induction frequency spans
                # multiple recovery delays
                if ind_freq == 50:
                    col_metrics['stp_initial_50hz'].append(initial)
                    if amps[1] != 0:
                        paired_pulse_ratio.append(amps[2] / amps[1])
            if all([k in pulses for k in [1,6,7,8]]):
                induction = (np.median([amps[6], amps[7], amps[8]]) - amps[1]) / amp_90p
                collect_induction.append(induction)
                if ind_freq == 50:
                    col_metrics['stp_induction_50hz'].append(induction)
            if delay is not None and all([k in pulses for k in range(1,13)]):
                r = [amps[i+8] - amps[i] for i in range(1,5)]
                recovery = np.median(r) / amp_90p
                collect_recovery.append(recovery)
                if delay == 250e-3:
                    col_metrics['stp_recovery_250ms'].append(recovery)
            if delay is not None and all([k in pulses for k in range(1,10)]):
                r = amps[9] - amps[1]
                recovery = np.median(r) / amp_90p
                collect_recovery_single.append(recovery)
                if delay == 250e-3:
                    col_metrics['stp_recovery_single_250ms'].append(recovery)

            # collect individual pulse amplitudes
            for i in range(1, 13):
                if i not in pulses:
                    break
                collect_pulse_amps[i-1].append(amps[i])
        
        stp_metrics = {
            'stp_initial': (np.median(collect_initial), np.std(collect_initial), len(collect_initial),) if len(collect_initial) > 1 else (float('nan'),)*3,
            'stp_induction': (np.median(collect_induction), np.std(collect_induction), len(collect_induction),) if len(collect_induction) > 1 else (float('nan'),)*3,
            'stp_recovery': (np.median(collect_recovery), np.std(collect_recovery), len(collect_recovery),) if len(collect_recovery) > 1 else (float('nan'),)*3,
            'stp_recovery_single': (np.median(collect_recovery_single), np.std(collect_recovery_single), len(collect_recovery_single),) if len(collect_recovery_single) > 1 else (float('nan'),)*3,
            'pulse_amplitudes': [
                (np.median(pulses), np.std(pulses), len(pulses),) if len(pulses) > 1 else float('nan')
                for i,pulses in enumerate(collect_pulse_amps)
            ]
        }
        all_metrics.append((meta, stp_metrics))
    
    # set one column to the full set of STP analysis
    setattr(dynamics, 'stp_all_stimuli', all_metrics)        
    
    # PPR is a bit out of place here, but we're including it since it's a popular metric used
    # in the literature.
    dynamics.paired_pulse_ratio_50hz = np.median(paired_pulse_ratio)
    
    # set individual columns for 50hz and 250ms
    for k,v in col_metrics.items():
        setattr(dynamics, k, np.median(v))
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
        return np.log((np.std(x)**2 - noise_std**2)**0.5 / abs(np.median(x)))

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
    amps_1_2 = ([], [])
    amps_2_4 = ([], [])
    amps_4_8 = ([], [])
    pulse_amps = [[] for i in range(9)]
    for key,recs in sorted_prs.items():
        clamp_mode, ind_freq, rec_delay = key
        if ind_freq != 50:
            continue
        for recording, pulses in recs.items():
            for i in range(1,9):
                if i not in pulses:
                    break
                if i == 2:
                    amps_1_2[0].append(getattr(pulses[i-1].PulseResponseFit, amp_field))
                    amps_1_2[1].append(getattr(pulses[i].PulseResponseFit, amp_field))
                if i == 4:
                    first = np.median([getattr(pulses[i].PulseResponseFit, amp_field) for i in range(1, 3)])
                    second = np.median([getattr(pulses[i].PulseResponseFit, amp_field) for i in range(3, 5)])
                    amps_2_4[0].append(np.median(first))
                    amps_2_4[1].append(np.median(second))
                if i == 8:
                    first = np.median([getattr(pulses[i].PulseResponseFit, amp_field) for i in range(1, 5)])
                    second = np.median([getattr(pulses[i].PulseResponseFit, amp_field) for i in range(5, 9)])
                    amps_4_8[0].append(np.median(first))
                    amps_4_8[1].append(np.median(second))
    
    if len(amps_1_2[0]) > 3:
        r,p = scipy.stats.pearsonr(amps_1_2[0], amps_1_2[1])
        dynamics.paired_event_correlation_1_2_r = r
        dynamics.paired_event_correlation_1_2_p = p
    if len(amps_2_4[0]) > 3:
        r,p = scipy.stats.pearsonr(amps_2_4[0], amps_2_4[1])
        dynamics.paired_event_correlation_2_4_r = r
        dynamics.paired_event_correlation_2_4_p = p
    if len(amps_4_8[0]) > 3:
        r,p = scipy.stats.pearsonr(amps_4_8[0], amps_4_8[1])
        dynamics.paired_event_correlation_4_8_r = r
        dynamics.paired_event_correlation_4_8_p = p



    return dynamics


def stim_sorted_pulse_amp(pair):
    qc_field = pair.synapse.synapse_type + '_qc_pass'

    q = db.query(
        db.PulseResponseFit.fit_amp,
        db.PulseResponseFit.dec_fit_reconv_amp,
        db.PulseResponseFit.baseline_dec_fit_reconv_amp,

        getattr(db.PulseResponse, qc_field).label('qc_pass'),
        db.StimPulse.pulse_number,
        db.MultiPatchProbe.induction_frequency,
        db.MultiPatchProbe.recovery_delay,
        db.SyncRec.ext_id.label('sync_rec_ext_id'),
    )
    q = q.join(db.PulseResponse, db.PulseResponseFit.pulse_response)
    q = q.join(db.Recording, db.PulseResponse.recording)
    q = q.join(db.SyncRec, db.Recording.sync_rec)
    q = q.join(db.PatchClampRecording, db.Recording.patch_clamp_recording)
    q = q.join(db.MultiPatchProbe, db.PatchClampRecording.multi_patch_probe)
    q = q.join(db.StimPulse, db.PulseResponse.stim_pulse)
    q = q.join(db.Experiment, db.SyncRec.experiment)
    q = q.join(db.Pair, db.PulseResponse.pair)
    q = q.join(db.Synapse, db.Pair.synapse)
    q = q.filter(db.Pair.id==pair.id)
    q = q.filter(db.PatchClampRecording.clamp_mode=='ic')
    q = q.order_by(db.SyncRec.ext_id, db.StimPulse.pulse_number)

    data = q.dataframe()

    qc_pass_data = data[data['qc_pass']]
    
    return qc_pass_data

def stp_all_stim_to_df(pairs, stp_df=None, pair_data=None):
    """
    Unpack dynamics.stp_all_stimuli column into dataframe for easier analysis

    Parameters
    ----------
    pairs: list
            Synphys Database pair items
    
    stp_df: Pandas Dataframe
            An exisiting dataframe to add this data to

    pair_data: dictionary
            pair metadata that you want in the resulting dataframe (ex. pre_cell_class or post_cell_class).
            Keys will be columns and values row data  

    Output
    -------
    stp_df: Pandas dataframe with stp data from all stimuli. Each row is a single induction frequency and
            recovery delay data point
    """
    metric = ['_median', '_std', '_n']
    stp_df = pd.Dataframe() if stp_df is None else stp_df
    for pair in pairs:
        if pair.dynamics is None:
            continue
        pair_data = {} if pair_data is None else pair_data
        pair_data.update({'pair_id': pair.id})
        stp_all_stim = pair.dynamics.stp_all_stimuli
        for stp in stp_all_stim:
            meta, data = stp
            pair_data.update({'clamp_mode': meta[0], 'ind_freq': meta[1], 'rec_delay': meta[2]})
            data2 = {}
            for k, value in data.items():
                if type(value) != list:
                    continue
                for m, v in zip(metric, value):
                    data2[k+m] = v
            pair_data.update(data2)
            stp_df = stp_df.append(pair_data, ignore_index=True)
    return stp_df