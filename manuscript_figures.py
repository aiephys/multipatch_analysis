import pickle
import numpy as np
from synaptic_dynamics import DynamicsAnalyzer
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.baseline import float_mode


def cache_response(expt, pre, post, cache_file, cache):
        key = (expt.uid, pre, post)
        if key in cache:
            responses = cache[key]
            amp_trials = len(responses['amp']['data'])
            amp = []
            if amp_trials != 0:
                for trial in range(amp_trials):
                    amp.append(Trace(data=responses['amp']['data'][trial], dt=responses['amp']['dt'][trial]))
            kinetics_trials = len(responses['kinetics']['data'])
            kinetics = []
            if kinetics_trials != 0:
                for trial in range(kinetics_trials):
                    kinetics.append(Trace(data=responses['kinetics']['data'][trial], dt=responses['kinetics']['dt'][trial]))
            return amp, kinetics

        amp, kinetics = get_response(expt, pre, post)
        cache[key] = {'amp': amp, 'kinetics': kinetics}

        data = pickle.dumps(cache)
        open(cache_file, 'wb').write(data)
        print "cached experiment %s" % key
        return amp, kinetics


def get_response(expt, pre, post):
    kinetics = {'data': [], 'dt': []}
    amp = {'data': [], 'dt': []}
    pulse = 0 # only pull first pulse
    analyzer = DynamicsAnalyzer(expt, pre, post, align_to='spike')
    pulse_response = analyzer.pulse_responses
    if len(pulse_response) == 0:
        print "No suitable data found for cell %d -> cell %d in expt %s" % (pre, post, expt.source_id)
        return amp, kinetics
    else:
        for i,stim_params in enumerate(pulse_response.keys()):
            resp = pulse_response[stim_params]
            ind_freq, rec_delay, holding = stim_params
            for trial in resp:
                r = trial[pulse]['response']
                if ind_freq <= 100:
                    amp['data'].append(r.data)
                    amp['dt'].append(r.dt)
                if ind_freq <= 20:
                    kinetics['data'].append(r.data)
                    kinetics['dt'].append(r.dt)
        return amp, kinetics

def get_amplitude(response_list):
    avg_trace = TraceList(response_list).mean()
    data = avg_trace.data
    dt = avg_trace.dt
    base = float_mode(data[:int(10e-3 / dt)])
    bsub_mean = avg_trace.copy(data=data - base)
    neg = data[int(13e-3/dt):].min() - base
    pos = data[int(13e-3/dt):].max() - base
    avg_amp = neg if abs(neg) > abs(pos) else pos
    amp_sign = '-' if avg_amp < 0 else '+'
    return bsub_mean, avg_amp, amp_sign

def align_to_rise(avg_trace, avg_amp):
    trace = avg_trace.data
    dt = avg_trace.dt
    base = trace[:int(10e-3/dt)]
    norm_amp = (trace - base)/avg_amp
    t = np.isclose(norm_amp, 1, atol=0.005)
    rise_times = np.where(t)
    if len(rise_times[0]) != 0:
        if rise_times[0].max() > 13e-3/dt:
            rise_time = rise_times[0][rise_times[0]>13e-3/dt][0]
            psp_align = trace.copy(data=trace.data[(rise_time - int(10e-3/dt)):])
            return psp_align
        else:
            return avg_trace
    else:
        return avg_trace