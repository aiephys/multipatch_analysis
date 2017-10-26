import pickle
import numpy as np
import pyqtgraph as pg
from synaptic_dynamics import DynamicsAnalyzer
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.baseline import float_mode
from neuroanalysis.event_detection import exp_deconvolve
from neuroanalysis.filter import bessel_filter
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison
app = pg.mkQApp()

def cache_response(expt, pre, post, cache_file, cache, type='pulse'):
        key = (expt.uid, pre, post)
        if key in cache:
            responses = cache[key]
            if type == 'pulse':
                response = format_responses(responses)
            else:
                response = responses
            return response

        responses = get_response(expt, pre, post, type=type)
        cache[key] = responses

        data = pickle.dumps(cache)
        open(cache_file, 'wb').write(data)
        print ("cached connection %s, %d -> %d" % (key[0], key[1], key[2]))
        if type == 'pulse':
            response = format_responses(responses)
        else:
            response = responses
        return response


def format_responses(responses):
    n_trials = len(responses['data'])
    response = {}
    if n_trials != 0:
        for trial in range(n_trials):
            stim_params = responses['stim_param'][trial]
            if stim_params not in response:
                response[stim_params] = []
            response[stim_params].append(Trace(data=responses['data'][trial], dt=responses['dt'][trial],
                                    stim_param=[responses['stim_param'][trial]]))
    return response

def get_response(expt, pre, post, type='pulse'):
    analyzer = DynamicsAnalyzer(expt, pre, post, method='deconv', align_to='spike')
    if type == 'pulse':
        pulse = 0  # only pull first pulse
        response = {'data': [], 'dt': [], 'stim_param': []}
        responses = analyzer.pulse_responses
        for i,stim_params in enumerate(responses.keys()):
            resp = responses[stim_params]
            for trial in resp:
                r = trial[pulse]['response']
                r.meta['stim_params'] = stim_params
                response['data'].append(r.data)
                response['dt'].append(r.dt)
                response['stim_param'].append(r.meta['stim_params'])
    elif type == 'train':
        responses = analyzer.train_responses
        pulse_offset = analyzer.pulse_offsets
        response = {'responses': responses, 'pulse_offsets': pulse_offset}
    else:
        "Must select either pulse responses or train responses"
    if len(responses) == 0:
        print "No suitable data found for cell %d -> cell %d in expt %s" % (pre, post, expt.source_id)
        return response

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
    bsub_mean = bsub(avg_trace)
    return bsub_mean

def bsub(trace):
    data = trace.data
    dt = trace.dt
    base = float_mode(data[:int(10e-3 / dt)])
    bsub_trace = trace.copy(data=data - base)
    return bsub_trace

def response_filter(response, freq_range=None, holding_range=None, pulse=False, train=None, delta_t=None):
    #plot = pg.plot()
    new_responses = []
    for stim_params, trials in response.items():
        ind_freq, rec_t, holding = stim_params
        holding = holding * 1e3
        rec_t = int(np.round(rec_t * 1e3, -1))
        if freq_range is not None and (ind_freq < freq_range[0] or ind_freq > freq_range[1]):
            continue
        if holding_range is not None and (holding > holding_range[0] or holding < holding_range[1]):
            continue
        if delta_t is not None and rec_t != delta_t:
            continue
        if pulse is True:
            for trial in trials:
                new_responses.append(trial)
        elif train is not None:
            for trial in trials[train].responses:
                new_responses.append(trial)
        else:
            new_responses.append(trials)
        # plot.plot(response[trial].time_values, response[trial].data)
        # app.processEvents()
    return new_responses

def feature_anova(feature, data):
    feature_list = [data[group][feature] for group in data.keys()]
    f, p = stats.f_oneway(feature_list[0], feature_list[1], feature_list[2], feature_list[3])
    print ('One-way ANOVA: %s' % feature)
    print ('=============')

    print ('F value: %.3f' % f)
    print ('P value: %.5f \n' % p)

def trace_plot(trace, color, plot=None, x_range=None):
    if plot is None:
        plot = pg.plot()
        plot.setLabels(left=('Vm', 'V'))
        plot.set_labels(bottom=('t', 's'))
        plot.setXRange(x_range)
    plot.plot(trace.time_values, trace.data, pen=color)
    app.processEvents()
    return plot

def train_amp(trace, pulse_offset, sign):
    deconv_trace = deconv_train(trace)
    pulses = np.array(pulse_offset)[0]
    ind_pulses = pulses[:8] + 10e-3
    rec_pulses = pulses[8:].copy()
    rec_pulses += 10e-3 - rec_pulses[0]
    ind = deconv_trace[0]
    rec = deconv_trace[1]

    amps = np.empty((len(ind), len(pulses)))
    for k, pulse_part, trace_part in [(0, ind_pulses, ind), (1, rec_pulses, rec)]:
        for i, n in enumerate(trace_part):
            dt = n.dt
            for j, pulse in enumerate(pulse_part):
                if k == 1:
                    j += 8
                start = int(pulse/dt)
                stop = start + int(4e-3/dt)
                chunk = n.data[start:stop]
                if sign == '+':
                    imx = np.argmax(chunk)
                else:
                    imx = np.argmin(chunk)
                mx = chunk[imx]
                amps[i, j] = mx

    return amps

def deconv_train(trace):
    deconv = [[], []]
    for i, k in enumerate(trace):
        for n in k:
            n_dec = bessel_filter(exp_deconvolve(n, 15e-3), 500)
            deconv[i].append(n_dec)

    return deconv