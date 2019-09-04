"""Quality control measurements

QC functions meant to ensure consistent filtering across different analyses
"""
import numpy as np
import pyqtgraph as pg
from neuroanalysis.util.data_test import DataTestCase


def recording_qc_pass(rec):
    """Applies a minimal set of QC criteria to a recording:

    * Must be a complete sweep (cannot contain large chunks of 0s)
    * Baseline RMS noise must be < 5mV or < 200 pA
    * Baseline current must be < 800 pA
    * For current clamp, baseline potential must be between -45 and -85 mV

    This is intended only to remove the most egregious data -- cells that are dead,
    sweeps that were interrupted before completion, etc. This is NOT intended to
    detect unhealthy cells, bad access resistance, etc.

    Parameters
    ----------
    rec : PatchClampRecording
        The PatchClampRecording instance to evaluate

    Returns
    -------
    qc_pass: bool
        Whether this recording passes qc_pass
    failures: list
        qc failures
    """
    failures = []

    if rec.baseline_current is None:
       failures.append('unknown baseline current')
    elif rec.baseline_current < -800e-12 or rec.baseline_current > 800e-12:
       failures.append('baseline current of %s is outside of bounds [-800pA, 800pA]' % pg.siFormat(rec.baseline_current, suffix='A'))
    
    if rec.clamp_mode == 'ic':
        if rec.baseline_potential is None:
            failures.append('baseline potential is None')
        elif rec.baseline_potential < -85e-3 or rec.baseline_potential > -45e-3:
            failures.append('baseline potential of %s is outside of bounds [-85mV, -45mV]' % pg.siFormat(rec.baseline_potential, suffix='V'))
        
        if rec.baseline_rms_noise is None:
            failures.append('no baseline_rms_noise for this recording')
        elif rec.baseline_rms_noise > 5e-3:
            failures.append('baseline rms noise of %s exceeds 5mV' % pg.siFormat(rec.baseline_rms_noise, suffix='V'))
        
    elif rec.clamp_mode == 'vc':
        if rec.baseline_rms_noise is None:
            failures.append('no baseline_rms_noise for this recording')
        elif rec.baseline_rms_noise > 200e-12:
           failures.append('baseline rms noise of %s exceeds 200pA' % pg.siFormat(rec.baseline_rms_noise, suffix='A'))
       
        
    data = rec['primary'].data
    if (data == 0).sum() > len(data) // 10:
        failures.append('data recording contains a significant chunk of zeros')

    qc_pass = len(failures)==0
    return qc_pass, failures


def pulse_response_qc_pass(post_rec, window, n_spikes, adjacent_pulses):
    """Apply QC criteria for pulse-response recordings:

    * Postsynaptic recording passes recording_qc_pass()
    * Presynaptic cell must have at least 1 spike in response to pulse
    * No other presynaptic pulses within 8ms on either side
    * Inhibitory response baseline potential must be between -45 and -60 mV
    * Excitatory response baseline potential must be between -45 and -80 mV
    * Overall stdev for postsynaptic recording must be < 1.5 mV or < 15 pA
    * Current clamp response must never exceed -40 mV
    
    These criteria are intended as minimal quality control when determining _whether_ a synaptic
    connection exists between two cells. 

    Parameters
    ----------
    post_rec : Recording
        The postsynaptic Recording instance
    window : list
        [start, stop] times indicating the region of the postsynaptic recording containing the pulse response
    n_spikes : int or None
        The number of presynaptic spikes evoked for this pulse response. If None, then this
        check is skipped (this is used for background data where we do not expect to have spikes).
    adjacent_pulses : list
        The times of any adjacent presynaptic stimulus pulses, relative to the spike of interest.
        This is used to ensure there is a minimum window of quiescence around the pulse to test, which
        excludes responses from very high frequency stimuli.

    Returns
    -------
    ex_qc_pass : bool
        Whether this pulse-response passes QC for detecting excitatory connections
    in_qc_pass : bool
        Whether this pulse-response passes QC for detecting inhibitory connections
    failures : dict
        QC failures for ex and in
    """
    
    failures = {'ex': [], 'in': []}

    # Require the postsynaptic recording to pass basic QC
    recording_pass_qc, recording_qc_failures = recording_qc_pass(post_rec)
    if recording_pass_qc is False:
        [failures[k].append('postsynaptic recording failed QC: %s' % ', and '.join(recording_qc_failures)) for k in failures.keys()]

    # require at least 1 presynaptic spike
    if n_spikes == 0:
        [failures[k].append('%d spikes detected in presynaptic recording' % n_spikes) for k in failures.keys()]

    # Check for noise in response window
    data = post_rec['primary'].time_slice(window[0], window[1])
    pre_pulse = data.time_slice(window[0], window[0] + 5e-3)
    base = pre_pulse.median()
    max_amp = np.abs((data - base).data).max()
    if post_rec.clamp_mode == 'ic':
        base_potential = base
        if pre_pulse.std() > 1.5e-3:
            [failures[k].append('STD of response window, %s, exceeds 1.5mV' % pg.siFormat(pre_pulse.std(), suffix='V')) for k in failures.keys()]
        if data.data.max() > -40e-3:
            [failures[k].append('Max in response window, %s, exceeds -40mV' % pg.siFormat(data.data.max(), suffix='V')) for k in failures.keys()]
        if max_amp > 10e-3:
            [failures[k].append('Max response amplitude, %s, exceeds 10mV' % pg.siFormat(max_amp, suffix='V')) for k in failures.keys()]
    elif post_rec.clamp_mode == 'vc':
        base_potential = post_rec['command'].time_slice(window[0], window[1]).median()
        if pre_pulse.std() > 15e-12:
            [failures[k].append('STD of response window, %s, exceeds 15pA' % pg.siFormat(pre_pulse.std(), suffix='A')) for k in failures.keys()]
        if max_amp > 500e-12:
            [failures[k].append('Max response amplitude, %s, exceeds 500pA' % pg.siFormat(max_amp, suffix='A')) for k in failures.keys()]
    else:
        raise TypeError('Unsupported clamp mode %s' % post_rec.clamp_mode)

    # Check timing of adjacent spikes
    if any([abs(t) < 8e-3 for t in adjacent_pulses]):
        [failures[k].append('Spikes detected within 8ms of the response window') for k in failures.keys()]

    # Check holding potential is appropriate for each sign
    ex_limits = [-85e-3, -45e-3]
    in_limits = [-60e-3, -45e-3]
    # check both baseline_potential (which is measured over all baseline regions in the recording)
    # and *base_potential*, which is just the median value over the IC pre_pulse window or VC command
    
    if not (ex_limits[0] < base_potential < ex_limits[1]): 
        failures['ex'].append('Response window baseline of %s is outside of bounds [-85mV, -45mV]' % pg.siFormat(base_potential, suffix='V'))
    if not (in_limits[0] < base_potential < in_limits[1]): 
        failures['in'].append('Response window baseline of %s is outside of bounds [-60mV, -45mV]' % pg.siFormat(base_potential, suffix='V'))
    
    base2 = post_rec.baseline_potential
    if base2 is None:
        failures['ex'].append('Unknown baseline potential for this recording')
        failures['in'].append('Unknown baseline potential for this recording')
    else:
        if not (ex_limits[0] < base2 < ex_limits[1]):
            failures['ex'].append('Recording baseline of %s is outside of bounds [-85mV, -45mV]' % pg.siFormat(base2, suffix='V'))
        if not (in_limits[0] < base2 < in_limits[1]):
            failures['in'].append('Recording baseline of %s is outside of bounds [-60mV, -45mV]' % pg.siFormat(base2, suffix='V'))
    
    
    ex_qc_pass = len(failures['ex'])==0 
    in_qc_pass = len(failures['in'])==0
    
    return ex_qc_pass, in_qc_pass, failures

def spike_qc(n_spikes, post_qc):
    """If there is not exactly 1 presynaptic spike, qc Fail spike and postsynaptic response
    """

    spike_qc_pass = n_spikes == 1
    trace_qc_pass = False if spike_qc_pass is False else post_qc

    return spike_qc_pass, trace_qc_pass

class PulseResponseQCTestCase(DataTestCase):
    def __init__(self):
        DataTestCase.__init__(self, pulse_response_qc_pass)

    @property
    def name(self):
        meta = self.meta
        return "%s_%s_%s_%0.3f" % (meta['expt_id'], meta['sweep_id'], meta['post_cell_id'], self.input_args['window'][0])