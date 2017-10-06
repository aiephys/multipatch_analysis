import pickle
import numpy as np
from synaptic_dynamics import DynamicsAnalyzer
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.baseline import float_mode
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison


def cache_response(expt, pre, post, cache_file, cache):
        key = (expt.uid, pre, post)
        if key in cache:
            responses = cache[key]
            n_trials = len(responses['data'])
            response = []
            if n_trials != 0:
                for trial in range(n_trials):
                    response.append(Trace(data=responses['data'][trial], dt=responses['dt'][trial],
                                        stim_param=[responses['stim_param'][trial]]))
            return response

        response = get_response(expt, pre, post)
        cache[key] = response

        data = pickle.dumps(cache)
        open(cache_file, 'wb').write(data)
        print ("cached connection %s, %d -> %d" % (key[0], key[1], key[2]))
        n_trials = len(response['data'])
        response2 = []
        if n_trials != 0:
            for trial in range(n_trials):
                response2.append(Trace(data=response['data'][trial], dt=response['dt'][trial],
                                  stim_param=[response['stim_param'][trial]]))
        return response2


def get_response(expt, pre, post):
    response = {'data': [], 'dt': [], 'stim_param': []}
    pulse = 0 # only pull first pulse
    analyzer = DynamicsAnalyzer(expt, pre, post, align_to='spike')
    pulse_response = analyzer.pulse_responses
    if len(pulse_response) == 0:
        print "No suitable data found for cell %d -> cell %d in expt %s" % (pre, post, expt.source_id)
        return response
    else:
        for i,stim_params in enumerate(pulse_response.keys()):
            resp = pulse_response[stim_params]
            for trial in resp:
                r = trial[pulse]['response']
                r.meta['stim_params'] = stim_params
                response['data'].append(r.data)
                response['dt'].append(r.dt)
                response['stim_param'].append(r.meta['stim_params'])
        return response

def get_amplitude(response_list):
    avg_trace = TraceList(response_list).mean()
    data = avg_trace.data
    dt = avg_trace.dt
    base = float_mode(data[:int(10e-3 / dt)])
    bsub_mean = avg_trace.copy(data=data - base)
    neg = bsub_mean.data[int(13e-3/dt):].min()
    pos = bsub_mean.data[int(13e-3/dt):].max()
    avg_amp = neg if abs(neg) > abs(pos) else pos
    amp_sign = '-' if avg_amp < 0 else '+'
    peak_ind = list(bsub_mean.data).index(avg_amp)
    peak_t = bsub_mean.time_values[peak_ind]
    return bsub_mean, avg_amp, amp_sign, peak_t

def response_filter(response, freq_range=None, holding_range=None):
    new_responses = []
    for trial in range(len(response)):
        ind_freq, _, holding = response[trial].meta['stim_param'][0]
        holding = holding * 1e3
        if freq_range is not None and (ind_freq < freq_range[0] or ind_freq > freq_range[1]):
            continue
        if holding is not None and (holding > holding_range[0] or holding < holding_range[1]):
            continue
        new_responses.append(response[trial])
    return new_responses

def feature_anova(feature, data):
    feature_list = [data[group][feature] for group in data.keys()]
    f, p = stats.f_oneway(feature_list[0], feature_list[1], feature_list[2], feature_list[3])
    print ('One-way ANOVA: %s' % feature)
    print ('=============')

    print ('F value: %.3f' % f)
    print ('P value: %.5f \n' % p)