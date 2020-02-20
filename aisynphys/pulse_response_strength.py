# coding: utf8
"""
Analyses that measure the strength of synaptic connections.

"""
from __future__ import print_function, division

import sys, multiprocessing, time, warnings

import numpy as np
import pyqtgraph as pg

from neuroanalysis.data import TSeries
from neuroanalysis import filter
from neuroanalysis.event_detection import exp_deconvolve, exp_reconvolve, exp_deconv_psp_params
from neuroanalysis.fitting import fit_psp, Psp, SearchFit, fit_scale_offset
from neuroanalysis.baseline import float_mode

from .database import default_db as db


def measure_response(pr):
    """Curve fit a single pulse response to measure its amplitude / kinetics.
    
    Uses the known latency and kinetics of the synapse to seed the fit.
    Optionally fit a baseline at the same time for noise measurement.
    
    Parameters
    ----------
    pr : PulseResponse
    """
    syn = pr.pair.synapse
    pcr = pr.recording.patch_clamp_recording
    if pcr.clamp_mode == 'ic':
        rise_time = syn.psp_rise_time
        decay_tau = syn.psp_decay_tau
    else:
        rise_time = syn.psc_rise_time
        decay_tau = syn.psc_decay_tau
                
    # make sure all parameters are available
    for v in [pr.stim_pulse.first_spike_time, syn.latency, rise_time, decay_tau]:
        if v is None or syn.latency is None or not np.isfinite(v):
            return None, None
    
    data = pr.get_tseries('post', align_to='spike')

    # decide whether/how to constrain the sign of the fit
    if syn.synapse_type == 'ex':
        sign = 1
    elif syn.synapse_type == 'in':
        if pcr.baseline_potential > -60e-3:
            sign = -1
        else:
            sign = 0
    else:
        sign = 0
    if pcr.clamp_mode == 'vc':
        sign = -sign

    # fit response region
    response_fit = fit_psp(data,
        search_window=syn.latency + np.array([-100e-6, 100e-6]), 
        clamp_mode=pcr.clamp_mode, 
        sign=sign,
        baseline_like_psp=True, 
        init_params={'rise_time': rise_time, 'decay_tau': decay_tau},
        refine=False,
        decay_tau_bounds=('fixed',),
        rise_time_bounds=('fixed',),        
    )
        
    # fit baseline region
    baseline = pr.get_tseries('baseline', align_to='spike')
    if baseline is None:
        baseline_fit = None
    else:
        baseline_fit = fit_psp(baseline,
            search_window=syn.latency + np.array([-100e-6, 100e-6]), 
            clamp_mode=pcr.clamp_mode, 
            sign=sign, 
            baseline_like_psp=True, 
            init_params={'rise_time': rise_time, 'decay_tau': decay_tau},
            refine=False,
        )

    return response_fit, baseline_fit


def measure_deconvolved_response(pr):
    """Use exponential deconvolution and a curve fit to estimate the amplitude of a synaptic response.

    Uses the known latency and kinetics of the synapse to constrain the fit.
    Optionally fit a baseline at the same time for noise measurement.
    
    Parameters
    ----------
    pr : PulseResponse
    """
    syn = pr.pair.synapse
    pcr = pr.recording.patch_clamp_recording
    if pcr.clamp_mode == 'ic':
        rise_time = syn.psp_rise_time
        decay_tau = syn.psp_decay_tau
        lowpass = 2000
    else:
        rise_time = syn.psc_rise_time
        decay_tau = syn.psc_decay_tau
        lowpass = 6000
    
    # make sure all parameters are available
    for v in [pr.stim_pulse.first_spike_time, syn.latency, rise_time, decay_tau]:
        if v is None or not np.isfinite(v):
            return None, None

    response_data = pr.get_tseries('post', align_to='spike')
    baseline_data = pr.get_tseries('baseline', align_to='spike')

    ret = []
    for data in (response_data, baseline_data):
        if data is None:
            ret.append(None)
            continue
            
        filtered = deconv_filter(data, None, tau=decay_tau, lowpass=lowpass, remove_artifacts=False, bsub=True)
        
        # chop down to the minimum we need to fit the deconvolved event.
        # there's a tradeoff here -- to much data and we risk incorporating nearby spontaneous events; too little
        # data and we get more noise in the fit to baseline
        filtered = filtered.time_slice(syn.latency-1e-3, syn.latency + rise_time + 1e-3)
        
        # Deconvolving a PSP-like shape yields a narrower PSP-like shape with lower rise power.
        # Guess the deconvolved time constants:
        dec_amp, dec_rise_time, dec_rise_power, dec_decay_tau = exp_deconv_psp_params(amp=1, rise_time=rise_time, decay_tau=decay_tau, rise_power=2)
        amp_ratio = 1 / dec_amp
        
        psp = Psp()

        # Need to measure amplitude of exp-deconvolved events; two methods to pick from here:
        # 1) Direct curve fitting using the expected deconvolved rise/decay time constants. This
        #    allows some wiggle room in latency, but produces a weird butterfly-shaped background noise distribution.
        # 2) Analytically calculate the scale/offset of a fixed template. Uses a fixed latency, but produces
        #    a nice, normal-looking background noise distribution.

        # Measure amplitude of deconvolved event by curve fitting:
        # with warnings.catch_warnings():
        #     warnings.simplefilter("ignore")
        #     max_amp = filtered.data.max() - filtered.data.min()
        #     fit = psp.fit(filtered.data, x=filtered.time_values, params={
        #         'xoffset': (response_rec.latency, response_rec.latency-0.2e-3, response_rec.latency+0.5e-3),
        #         'yoffset': (0, 'fixed'),
        #         'amp': (0, -max_amp, max_amp),
        #         'rise_time': (dec_rise_time, 'fixed'),
        #         'decay_tau': (dec_decay_tau, 'fixed'),
        #         'rise_power': (dec_rise_power, 'fixed'),
        #     })
        # reconvolved_amp = fit.best_values['amp'] * amp_ratio
        
        # fit = {
        #     'xoffset': fit.best_values['xoffset'],
        #     'yoffset': fit.best_values['yoffset'],
        #     'amp': fit.best_values['amp'],
        #     'rise_time': dec_rise_time,
        #     'decay_tau': dec_decay_tau,
        #     'rise_power': dec_rise_power,
        #     'reconvolved_amp': reconvolved_amp,
        # }

        # Measure amplitude of deconvolved events by direct template match
        template = psp.eval(
            x=filtered.time_values, 
            xoffset=syn.latency,
            yoffset=0,
            amp=1,
            rise_time=dec_rise_time,
            decay_tau=dec_decay_tau,
            rise_power=dec_rise_power,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            scale, offset = fit_scale_offset(filtered.data, template)

        # calculate amplitude of reconvolved event -- tis is our best guess as to the
        # actual event amplitude
        reconvolved_amp = scale * amp_ratio
        
        fit = {
            'xoffset': syn.latency,
            'yoffset': offset,
            'amp': scale,
            'rise_time': dec_rise_time,
            'decay_tau': dec_decay_tau,
            'rise_power': 1,
            'reconvolved_amp': reconvolved_amp,
        }
        
        ret.append(fit)
        
    return ret


def measure_peak(trace, sign, spike_time, pulse_times, spike_delay=1e-3, response_window=4e-3):
    # Start measuring response after the pulse has finished, and no earlier than 1 ms after spike onset
    # response_start = max(spike_time + spike_delay, pulse_times[1])

    # Start measuring after spike and hope that the pulse offset doesn't get in the way
    # (if we wait for the pulse to end, then we miss too many fast rise / short latency events)
    response_start = spike_time + spike_delay
    response_stop = response_start + response_window

    # measure baseline from beginning of data until 50µs before pulse onset
    baseline_start = trace.t0
    baseline_stop = pulse_times[0] - 50e-6

    baseline = float_mode(trace.time_slice(baseline_start, baseline_stop).data)
    response = trace.time_slice(response_start, response_stop)
    if (response.t_end - response.t0) < 0.8 * response_window:
        # response window is too short; don't attempt to make a measurement.
        # (this only happens when the spike is very late and the next pulse is very soon)
        return None, None

    if sign == '+':
        i = np.argmax(response.data)
    else:
        i = np.argmin(response.data)
    peak = response.data[i]
    latency = response.time_values[i] - spike_time
    return peak - baseline, latency


def measure_sum(trace, sign, baseline=(0e-3, 9e-3), response=(12e-3, 17e-3)):
    baseline = trace.time_slice(*baseline).data.sum()
    peak = trace.time_slice(*response).data.sum()
    return peak - baseline

        
def deconv_filter(trace, pulse_times, tau=15e-3, lowpass=24000., lpf=True, remove_artifacts=False, bsub=True):
    if tau is not None:
        dec = exp_deconvolve(trace, tau)
    else:
        dec = trace

    if remove_artifacts:
        # after deconvolution, the pulse causes two sharp artifacts; these
        # must be removed before LPF
        cleaned = remove_crosstalk_artifacts(dec, pulse_times)
    else:
        cleaned = dec

    if bsub:
        baseline = np.median(cleaned.time_slice(cleaned.t0+5e-3, cleaned.t0+10e-3).data)
        b_subbed = cleaned - baseline
    else:
        b_subbed = cleaned

    if lpf:
        return filter.bessel_filter(b_subbed, lowpass)
    else:
        return b_subbed


def remove_crosstalk_artifacts(data, pulse_times):
    dt = data.dt
    r = [-50e-6, 250e-6]
    edges = [(int((t+r[0])/dt), int((t+r[1])/dt)) for t in pulse_times]
    # If window is too shortm then it becomes seneitive to sample noise.
    # If window is too long, then it becomes sensitive to slower signals (like the AP following pulse onset)
    return filter.remove_artifacts(data, edges, window=100e-6)


def analyze_response_strength(pr, source, remove_artifacts=False, deconvolve=True, lpf=True, bsub=True, lowpass=1000):
    """Perform a standardized strength analysis on a record selected by response_query.

    1. Determine timing of presynaptic stimulus pulse edges and spike
    2. Measure peak deflection on raw trace
    3. Apply deconvolution / artifact removal / lpf
    4. Measure peak deflection on deconvolved trace
    """
    if source == 'pulse_response':
        data = pr.get_tseries('post', align_to='pulse')
    elif source == 'baseline':
        data = pr.get_tseries('baseline', align_to='pulse')
    else:
        raise ValueError("Invalid source %s" % source)
    
    if data is None:
        return None
            
    pulse_times = data.meta['pulse_start'], data.meta['pulse_stop']
    spike_time = data.meta['spike_time']
    if spike_time is None:
        # these pulses failed QC, but we analyze them anyway to make all data visible
        spike_time = 1e-3

    results = {}

    results['raw_trace'] = data
    results['pulse_times'] = pulse_times
    results['spike_time'] = spike_time

    # Measure crosstalk from pulse onset
    p1 = data.time_slice(pulse_times[0]-200e-6, pulse_times[0]).median()
    p2 = data.time_slice(pulse_times[0], pulse_times[0]+200e-6).median()
    results['crosstalk'] = p2 - p1

    # crosstalk artifacts in VC are removed before deconvolution
    pcr = pr.recording.patch_clamp_recording
    if pcr.clamp_mode == 'vc' and remove_artifacts is True:
        data = remove_crosstalk_artifacts(data, pulse_times)
        remove_artifacts = False

    # Measure deflection on raw data
    results['pos_amp'], _ = measure_peak(data, '+', spike_time, pulse_times)
    results['neg_amp'], _ = measure_peak(data, '-', spike_time, pulse_times)
    
    # Deconvolution / artifact removal / filtering
    if deconvolve:
        tau = 15e-3 if pcr.clamp_mode == 'ic' else 5e-3
    else:
        tau = None
    dec_data = deconv_filter(data, pulse_times, tau=tau, lpf=lpf, remove_artifacts=remove_artifacts, bsub=bsub, lowpass=lowpass)
    results['dec_trace'] = dec_data

    # Measure deflection on deconvolved data
    results['pos_dec_amp'], results['pos_dec_latency'] = measure_peak(dec_data, '+', spike_time, pulse_times)
    results['neg_dec_amp'], results['neg_dec_latency'] = measure_peak(dec_data, '-', spike_time, pulse_times)
    
    return results
