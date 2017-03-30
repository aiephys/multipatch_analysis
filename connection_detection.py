def detect_pulses(trace):
    sdiff = np.diff(trace.data)
    on_times = np.argwhere(sdiff > 0)[1:, 0]  # 1: skips test pulse
    off_times = np.argwhere(sdiff < 0)[1:, 0]
    return on_inds, off_inds


def detect_patch_evoked_spikes(pre_rec):
    """Given presynaptic Recording, detect action potentials
    evoked by current injection or unclamped spikes evoked by a voltage pulse.
    """
    pre_trace = pre_rec['primary']
    stim = pre_rec['command'].data

    # Detect pulse times
    on_inds, off_inds = detect_pulses(stim)

    # detect spike times
    spike_info = []
    for on, off in zip(on_inds, off_inds):
        spike = detect_evoked_spike(pre_rec, [on, off])
        spike_info.append(spike)
    return spike_info
    


def detect_connections(expt, max_freq):
    """
    loop over all sweeps (presynaptic)
        ignore sweeps with high induction frequency
        detect presynaptic spikes
        loop over all other sweeps in the same recording (postsynaptic)
            ignore sweeps that fail QC
            collect data following each stimulus
            collect baseline noise data
                (use noise.std() / sqrt(n) to estimate noise for averaged traces)
            collect mode / holding
    Average together all collected sweep chunks from the same cell, and of the same holding and clamp mode
        all pulses
        pulses 1, 2, 9, 10
        pulses 7, 8, 11, 12
    Fit result to psp
    Compare fit amplitude to noise floor
    Z-score fit NRMSE against fits from pre-stimulus data
        correct for multiple comparisons?
    Additional metric to detect low release probability connections?
    """
    expt_data = expt.data
    
    # loop over all sweeps (presynaptic)
    for srec in expt_data.contents:
        for pre_rec in srec.recordings:
            # todo: ignore sweeps with high induction frequency
            
            # detect presynaptic spikes
            spikes = detect_patch_evoked_spikes(pre_rec)
            
            # loop over all other sweeps in the same recording (postsynaptic)
            for post_rec in srec.recordings:
                if post_rec is pre_rec:
                    continue

                