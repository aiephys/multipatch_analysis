"""Quality control measurements

QC functions meant to ensure consistent filtering across different analyses
"""
import numpy as np


def recording_qc_pass(rec):
    """Applies a minimal set of QC criteria to a recording:

    * Must be a complete sweep (cannot contain large chunks of 0s)
    * Baseline RMS noise must be < 5mV or < 200 pA
    * Baseline current must be < 800 pA
    * For current clamp, baseline potential must be between -45 and -85 mV

    This is intended only to remove the most eggregious data -- cells that are dead,
    sweeps that were interrupted before completion, etc. This is NOT intended to
    detect unhealthy cells, bad access resistance, etc.

    Parameters
    ----------
    rec : PatchClampRecording
        The PatchClampRecording instance to evaluate
    """
    if rec.baseline_current < -800e-12 or rec.baseline_current > 800e-12:
        return False
    if rec.clamp_mode == 'ic':
        if rec.baseline_potential < -85e-3 or rec.baseline_potential > -45e-3:
            return False
        if rec.baseline_rms_noise > 5e-3:
            return False
    elif rec.clamp_mode == 'vc':
        if rec.baseline_rms_noise > 200e-12:
            return False
        
    data = rec['primary'].data
    if (data == 0).sum() > len(data) // 10:
        return False

    return True


def pulse_response_qc_pass(sign, post_rec, window, n_spikes):
    """Apply QC criteria for pulse-response recordings:

    * Postsynaptic recording passes recording_qc_pass()
    * Presynaptic cell must have at least 1 spike in response to pulse
    * Inhibitory response baseline potential must be between -45 and -60 mV
    * Excitatory response baseline potential must be between -45 and -80 mV
    * Overall stdev for postsynaptic recording must be < 1.5 mV or < 15 pA
    * Current clamp response must never exceed -40 mV

    Parameters
    ----------
    sign : -1 or +1
        Indicates whether to apply QC criteria for inhibitory (-1) or excitatory (+1) tests
    post_rec : Recording
        The postsynaptic Recording instance
    window : list
        [start, stop] indices indicating the region of the postsynaptic recording containing the pulse response
    n_spikes : int or None
        The number of presynaptic spikes evoked for this pulse response. If None, then this
        check is skipped (this is used for background data where we do not expect to have spikes).
    """
    if recording_qc_pass(post_rec) is False:
        return False

    if n_spikes == 0:
        return False
    
    if post_rec.clamp_mode == 'ic':
        data = post_rec['primary'][window[0]:window[1]]
        base = data.median()
        if data.std() > 1.5e-3:
            return False
        if data.data.max() > -40e-3:
            return False
    elif post_rec.clamp_mode == 'vc':
        data = post_rec['primary'][window[0]:window[1]]
        base = post_rec['command'][window[0]:window[1]].median()
        if data.std() > 15e-12:
            return False
    else:
        raise TypeError('Unsupported clamp mode %s' % post_rec.clamp_mode)
    
    if sign == 1:
        bmin, bmax = [-85e-3, -45e-3]
    elif sign == -1:
        bmin, bmax = [-60e-3, -45e-3]
    else:
        raise ValueError("sign must be -1 or +1")       

    # check both baseline_potential (which is measured over all baseline regions in the recording)
    # and *base*, which is just the median value over the response window
    base2 = post_rec.baseline_potential
    if not ((bmin < base < bmax) and (bmin < base2 < bmax)):
        return False

    return True
