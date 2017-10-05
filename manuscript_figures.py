import pickle
import numpy as np
from synaptic_dynamics import DynamicsAnalyzer
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.baseline import float_mode
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison


def cache_response(expt, pre, post, cache_file, cache, type='pulse'):
        key = (expt.uid, pre, post)
        if key in cache:
            responses = cache[key]
            response = format_responses(responses, type)
            return response

        responses = get_response(expt, pre, post, type=type)
        cache[key] = responses

        data = pickle.dumps(cache)
        open(cache_file, 'wb').write(data)
        print ("cached connection %s, %d -> %d" % (key[0], key[1], key[2]))
        response = format_responses(responses, type)
        return response


def format_responses(responses, type):
    if type == 'pulse':
        n_trials = len(responses['data'])
        response = []
        if n_trials != 0:
            for trial in range(n_trials):
                response.append(Trace(data=responses['data'][trial], dt=responses['dt'][trial],
                                      stim_param=[responses['stim_param'][trial]]))
    elif type == 'train':
        n_trials = len(responses['data'][0])
        response = [[], []]
        if n_trials != 0:
            for block in range(len(responses['data'])):
                for trial in range(n_trials):
                    response[block].append(Trace(data=responses['data'][block][trial], dt=responses['dt'][block][trial],
                                                 stim_param=responses['stim_param'][block][trial]))
    return response

def get_response(expt, pre, post, type='pulse'):
    analyzer = DynamicsAnalyzer(expt, pre, post, align_to='spike')
    if type == 'pulse':
        pulse = 0  # only pull first pulse
        response = {'data': [], 'dt': [], 'stim_param': []}
        responses = analyzer.pulse_responses
    elif type == 'train':
        response = {'data': [[], []], 'dt': [[], []], 'stim_param': [[], []]} # {'data': [[ind], [rec]], ...}
        responses = analyzer.train_responses
    else:
        "Must select either pulse responses or train responses"
    if len(responses) == 0:
        print "No suitable data found for cell %d -> cell %d in expt %s" % (pre, post, expt.source_id)
        return response
    else:
        for i,stim_params in enumerate(responses.keys()):
            resp = responses[stim_params]
            if type == 'pulse':
                for trial in resp:
                    r = trial[pulse]['response']
                    r.meta['stim_params'] = stim_params
                    response['data'].append(r.data)
                    response['dt'].append(r.dt)
                    response['stim_param'].append(r.meta['stim_params'])
            elif type == 'train':
                for i, block in enumerate(resp):
                    for trial in block.responses:
                        trial.meta['stim_params'] = stim_params
                        response['data'][i].append(trial.data)
                        response['dt'][i].append(trial.dt)
                        response['stim_param'][i].append(trial.meta['stim_params'])

        return response

def get_amplitude(response_list):
    bsub_mean = trace_avg(response_list)
    dt = bsub_mean.dt
    neg = bsub_mean.data[int(13e-3/dt):].min()
    pos = bsub_mean.data[int(13e-3/dt):].max()
    avg_amp = neg if abs(neg) > abs(pos) else pos
    amp_sign = '-' if avg_amp < 0 else '+'
    peak_ind = list(bsub_mean.data).index(avg_amp)
    peak_t = bsub_mean.time_values[peak_ind]
    return bsub_mean, avg_amp, amp_sign, peak_t

def trace_avg(response_list):
    avg_trace = TraceList(response_list).mean()
    data = avg_trace.data
    dt = avg_trace.dt
    base = float_mode(data[:int(10e-3 / dt)])
    bsub_mean = avg_trace.copy(data=data - base)
    return bsub_mean

def response_filter(response, freq_range=None, holding_range=None, delta_t=None):
    new_responses = []
    for trial in range(len(response)):
        if type(response[trial].meta['stim_param']) is list:
            response[trial].meta['stim_param'] = response[trial].meta['stim_param'][0]
        ind_freq, rec_t, holding = response[trial].meta['stim_param']
        holding = holding * 1e3
        rec_t = int(np.round(rec_t * 1e3, -1))
        if freq_range is not None and (ind_freq < freq_range[0] or ind_freq > freq_range[1]):
            continue
        if holding_range is not None and (holding > holding_range[0] or holding < holding_range[1]):
            continue
        if delta_t is not None and rec_t != delta_t:
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