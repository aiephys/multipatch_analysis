import os
import pickle
import numpy as np
import pyqtgraph as pg
import seaborn as sns
import time
import sys
import datetime
import re
from aisynphys.synaptic_dynamics import DynamicsAnalyzer
from neuroanalysis.data import TSeries, TSeriesList
from neuroanalysis.baseline import float_mode
from neuroanalysis.event_detection import exp_deconvolve
from neuroanalysis.filter import bessel_filter
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.spike_detection import detect_ic_evoked_spike
from scipy import stats
from aisynphys.constants import EXCITATORY_CRE_TYPES, INHIBITORY_CRE_TYPES
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison
app = pg.mkQApp()

colors_human = [(247, 118, 118),
                (246, 197, 97),
                (100, 202, 103),
                (107, 155, 250),
                (162, 62, 247)]

colors_mouse = [(249, 144, 92),
                (100, 202, 103),
                (81, 221, 209),
                (45, 77, 247),
                (162, 62, 247)]

def arg_to_date(arg):
    if arg is None:
        return None
    parts = re.split('\D+', arg)
    return datetime.date(*map(int, parts))

def write_cache(cache, cache_file):
    print("writing cache to disk...")
    pickle.dump(cache, open(cache_file + '.new', 'wb'))
    if os.path.exists(cache_file):
        os.remove(cache_file)
    os.rename(cache_file + '.new', cache_file)
    print("Done!")

def load_cache(cache_file):
    if os.path.exists(cache_file):
        try:
            result_cache = pickle.load(open(cache_file, 'rb'))
            print ('loaded cache')
        except:
            result_cache = {}
            sys.excepthook(*sys.exc_info())
            print ('Error loading cache')
    else:
        result_cache = {}

    return result_cache

def cache_response(expt, pre, post, cache, analysis_type='pulse'):
        key = (expt.uid, pre, post)
        if key in cache:
            response = cache[key]
            # if analysis_type == 'pulse':
            #     response = format_responses(responses)
            # else:
            #     response = responses
            cache_change = 0
            return response, cache_change

        response = get_response(expt, pre, post, analysis_type=analysis_type)
        cache[key] = response
        cache_change = 1
        print ("cached connection %s, %d -> %d" % (key[0], key[1], key[2]))
        # if analysis_type == 'pulse':
        #     response = format_responses(responses)
        # else:
        #     response = responses
        return response, cache_change


def format_responses(responses):
    n_trials = len(responses['data'])
    response = {}
    if n_trials != 0:
        for trial in range(n_trials):
            stim_params = responses['stim_param'][trial]
            if stim_params not in response:
                response[stim_params] = []
            response[stim_params].append(TSeries(data=responses['data'][trial], dt=responses['dt'][trial],
                                    stim_param=[responses['stim_param'][trial]]))
    return response

def get_response(expt, pre, post, analysis_type='pulse'):
    analyzer = DynamicsAnalyzer(expt, pre, post, method='deconv', align_to='spike')
    if analysis_type == 'pulse':
        response = analyzer.pulse_responses
    elif analysis_type == 'train':
        responses = analyzer.train_responses
        pulse_offset = analyzer.pulse_offsets
        response = {'responses': responses, 'pulse_offsets': pulse_offset}
    else:
        "Must select either pulse responses or train responses"
    if len(response) == 0:
        print "No suitable data found for cell %d -> cell %d in expt %s" % (pre, post, expt.source_id)
        return response, None
    artifact = analyzer.cross_talk()
    return response, artifact

def get_amplitude(response_list):
    """
    Parameters
    ----------
    response_list : list of neuroanalysis.data.TSeriesView objects
        neuroanalysis.data.TSeriesView object contains waveform data. 
    """
    
    if len(response_list) == 1:
        bsub_mean = bsub(response_list[0])
    else:
        bsub_mean = trace_avg(response_list)
    dt = bsub_mean.dt
    neg = bsub_mean.data[int(13e-3/dt):].min()
    pos = bsub_mean.data[int(13e-3/dt):].max()
    avg_amp = neg if abs(neg) > abs(pos) else pos
    amp_sign = '-' if avg_amp < 0 else '+'
    peak_ind = list(bsub_mean.data).index(avg_amp)
    peak_t = bsub_mean.time_values[peak_ind]
    return bsub_mean, avg_amp, amp_sign, peak_t

def fail_rate(response_list, sign, peak_t):
    amps = []
    for response in response_list:
        bsub_mean = bsub(response)
        dt = bsub_mean.dt
        if sign == '+':
            amp = bsub_mean.data[int(13e-3/dt):].max()
        else:
            amp = bsub_mean.data[int(13e-3 / dt): int((peak_t + 3e-3) / dt)].max()
        amps.append(amp)
    return amps


def trace_avg(response_list):
# doc string commented out to discourage code reuse given the change of values of t0
#    """
#    Parameters
#    ----------
#    response_list : list of neuroanalysis.data.TSeriesView objects
#        neuroanalysis.data.TSeriesView object contains waveform data. 
#        
#    Returns
#    -------
#    bsub_mean : neuroanalysis.data.TSeries object
#        averages and baseline subtracts the ephys waveform data in the 
#        input response_list TSeriesView objects and replaces the .t0 value with 0. 
#    
#    """
    for trace in response_list: 
        trace.t0 = 0  #align traces for the use of TSeriesList().mean() funtion
    avg_trace = TSeriesList(response_list).mean() #returns the average of the wave form in a of a neuroanalysis.data.TSeries object 
    bsub_mean = bsub(avg_trace) #returns a copy of avg_trace but replaces the ephys waveform in .data with the base_line subtracted wave_form
    
    return bsub_mean

def bsub(trace):
    """Returns a copy of the neuroanalysis.data.TSeries object 
    where the ephys data waveform is replaced with a baseline 
    subtracted ephys data waveform.  
    
    Parameters
    ----------
    trace : neuroanalysis.data.TSeries object  
        
    Returns
    -------
    bsub_trace : neuroanalysis.data.TSeries object
       Ephys data waveform is replaced with a baseline subtracted ephys data waveform
    """
    data = trace.data # actual numpy array of time series ephys waveform
    dt = trace.dt # time step of the data
    base = float_mode(data[:int(10e-3 / dt)]) # baseline value for trace 
    bsub_trace = trace.copy(data=data - base) # new neuroanalysis.data.TSeries object for baseline subtracted data
    return bsub_trace

def response_filter(response, freq_range=None, holding_range=None, pulse=False, train=None, delta_t=None):
    new_responses = []
    holding_pass = []
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
                new_responses.append(trial[0]['response'])
            holding_pass.append(holding)
        elif train is not None:
            for trial in trials[train].responses:
                new_responses.append(trial)
        else:
            new_responses.append(trials)

    return new_responses

def feature_anova(feature, data):
    feature_list = [(key, group[feature])for key, group in data.items()]
    f, p = stats.f_oneway(feature_list[0][1], feature_list[1][1], feature_list[2][1], feature_list[3][1])
    print ('One-way ANOVA: %s' % feature)
    print ('=============')
    for i in feature_list:
        print ('%s: %.3f +- %.3f' % (i[0], np.mean(i[1])*1e3, stats.sem(i[1])*1e3))
    print ('F value: %.3f' % f)
    print ('P value: %.5f \n' % p)
    return feature_list

def feature_kw(feature, data):
    feature_list = [(key, group[feature]) for key, group in data.items()]
    h, p = stats.kruskal(feature_list[0][1], feature_list[1][1], feature_list[2][1], feature_list[3][1])

    print ('Kruskal-Wallace: %s' % feature)
    print ('=============')
    for i in feature_list:
        print ('%s: %.3f +- %.3f' % (i[0], np.median(i[1]) * 1e3, np.std(i[1]) * 1e3))
    print ('H value: %.3f' % h)
    print ('P value: %.5f \n' % p)
    return feature_list

def trace_plot(trace, color, plot=None, x_range=None, name=None):
    if plot is None:
        plot = pg.plot()
    plot.setLabels(left=('Vm', 'V'))
    plot.setLabels(bottom=('t', 's'))
    if x_range is not None:
        plot.setXRange(x_range[0], x_range[1])
    plot.plot(trace.time_values, trace.data, pen=color, name=name)
    app.processEvents()
    return plot

def train_amp(trace, pulse_offset, sign):
    deconv_trace = deconv_train(trace[:2])
    pulses = np.array(pulse_offset)[0]
    ind_pulses = pulses[:8] + 13e-3
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

def induction_summary(train_response, freqs, holding, thresh=5, ind_dict=None, offset_dict=None, uid=None):
    if ind_dict is None:
        ind_dict = {}
        offset_dict = {}
    for f, freq in enumerate(freqs):
        induction_traces = {}
        induction_traces['responses'] = response_filter(train_response['responses'], freq_range=[freq, freq],
                                                        holding_range=holding, train=0)
        induction_traces['pulse_offsets'] = response_filter(train_response['pulse_offsets'], freq_range=[freq, freq])
        ind_rec_traces = response_filter(train_response['responses'], freq_range=[freq, freq], holding_range=holding,
                                         train=1, delta_t=250)
        if len(induction_traces['responses']) >= thresh and len(ind_rec_traces) >= thresh:
            induction_avg = trace_avg(induction_traces['responses'])
            ind_rec_avg = trace_avg(ind_rec_traces)
            ind_rec_avg.t0 = induction_avg.time_values[-1] + 0.1
            if freq not in ind_dict.keys():
                ind_dict[freq] = [[], [], []]
            ind_dict[freq][0].append(induction_avg)
            ind_dict[freq][1].append(ind_rec_avg)
            ind_dict[freq][2].append(uid)
            offset_dict[freq] = induction_traces['pulse_offsets']

    return ind_dict, offset_dict

def recovery_summary(train_response, rec_t, holding, thresh=5, rec_dict=None, offset_dict=None, uid=None):
    if rec_dict is None:
        rec_dict = {}
        offset_dict = {}
    rec_ind_traces = response_filter(train_response['responses'], freq_range=[50, 50], holding_range=holding, train=0)
    for t, delta in enumerate(rec_t):
        recovery_traces = {}
        recovery_traces['responses'] = response_filter(train_response['responses'], freq_range=[50, 50],
                                                       holding_range=holding, train=1, delta_t=delta)
        #rec_ind_traces = response_filter(train_response['responses'], freq_range=[50, 50],
        #                                               holding_range=holding, train=0, delta_t=delta)
        recovery_traces['pulse_offsets'] = response_filter(train_response['pulse_offsets'], freq_range=[50, 50],
                                                           delta_t=delta)
        if len(recovery_traces['responses']) >= thresh:
            recovery_avg = trace_avg(recovery_traces['responses'])
            rec_ind_avg = trace_avg(rec_ind_traces)
            recovery_avg.t0 = (rec_ind_avg.time_values[-1]) + 0.1
            if delta not in rec_dict.keys():
                rec_dict[delta] = [[], [], []]
            rec_dict[delta][0].append(rec_ind_avg)
            rec_dict[delta][1].append(recovery_avg)
            rec_dict[delta][2].append(uid)
            offset_dict[delta] = recovery_traces['pulse_offsets']

    return rec_dict, offset_dict

def pulse_qc(responses, baseline=None, pulse=None, plot=None):
    """
    Parameters
    ----------
    responses : list of neuroanalysis.data.TSeriesView objects
        neuroanalysis.data.TSeriesView object contains waveform data. 
    base_line : float
        Factor by which to multiply the standard deviation of the baseline current.
    pulse : float
        Factor by which to multiply the standard deviation of the current during a pulse.
        Currently not in use.
    plot : pyqtgraph.PlotItem
        If not None, plot the data on the referenced pyqtgraph object.

    Returns
    ----------
    qc_pass : 
    """
    qc_pass = []
    avg, amp, _, peak_t = get_amplitude(responses)
    pulse_win = int((peak_t + 1e-3)/avg.dt)
    pulse_std = np.std(avg.data[pulse_win:])
    base_win = int(10e-3/avg.dt)
    base_std = np.std(avg.data[:base_win])
    for response in responses:
        response = bsub(response)
        data = response.data
        if np.mean(data[:base_win]) > (baseline * base_std):
            if plot is not None:
                plot.plot(response.time_values, response.data, pen='r')
        #deprecate? or make standard?
        # elif np.mean(data[pulse_win:]) > (pulse * pulse_std) and plot is not None:
        #    if plot is not None:
        #         plot.plot(response.time_values, response.data, pen='b')
        else:
            if plot is not None:
                plot.plot(response.time_values, response.data)
            qc_pass.append(response)
    if len(qc_pass) > 0:
        qc_trace = trace_avg(qc_pass)
        if plot is not None:
            plot.addLegend()
            plot.plot(qc_trace.time_values, qc_trace.data, pen={'color':'k', 'width':2}, name=('%d'% len(qc_pass)))
    return qc_pass

def train_qc(responses, offset, amp=None, sign=None, plot=None):
    qc_pass = [[],[], []]
    amps = train_amp(responses[:2], offset, sign=sign)
    for n in range(amps.shape[0]):
        if plot is not None:
            plot.plot(responses[0][n].time_values, responses[0][n].data, pen='r')
            plot.plot(responses[1][n].time_values, responses[1][n].data, pen='r')
        if sign == '+' and np.any(amps[n, :] < 0):
            continue
        if sign == '-' and np.any(amps[n, :] > 0):
            continue
        if abs(np.mean(amps[n, :])) < amp:
            continue
        qc_pass[0].append(responses[0][n])
        qc_pass[1].append(responses[1][n])
        qc_pass[2].append(responses[2][n])
        if plot is not None:
            plot.plot(responses[0][n].time_values, responses[0][n].data, pen=[0, 0, 0, 100])
            plot.plot(responses[1][n].time_values, responses[1][n].data, pen=[0, 0, 0, 100])
            app.processEvents()
            time.sleep(1)
    return qc_pass

def subplots(name=None, row=None):
    p1 = name.addPlot(row=row, col=0)
    p2 = name.addPlot(row=row, col=1)
    p3 = name.addPlot(row=row, col=2)
    p4 = name.addPlot(row=row, col=3)
    p5 = name.addPlot(row=row, col=4)
    if row == 0:
        p1.setTitle('First Pulse Response, Example Connection')
        p2.setTitle('First Pulse Response, All Connections')
        p3.setTitle('50 Hz Train Response')
        p4.setTitle('Induction Amplitudes')
        p5.setTitle('Recovery Amplitudes')
    return p1, p2, p3, p4, p5


def summary_plot_pulse(feature_list, labels, titles, i, median=False, grand_trace=None, plot=None, color=None, name=None):
    """
    Plots features of single-pulse responses such as amplitude, latency, etc. for group analysis. Can be used for one
    group by ideal for comparing across many groups in the feature_list

    Parameters
    ----------
    feature_list : list of lists of floats
        single-pulse features such as amplitude. Can be multiple features each a list themselves
    labels : list of pyqtgraph.LabelItem
        axis labels, must be a list of same length as feature_list
    titles : list of strings
        plot title, must be a list of same length as feature_list
    i : integer
        iterator to place groups along x-axis
    median : boolean
        to calculate median (True) vs mean (False), default is False
    grand_trace : neuroanalysis.data.TSeriesView object
        option to plot response trace alongside scatter plot, default is None
    plot : pyqtgraph.PlotItem
        If not None, plot the data on the referenced pyqtgraph object.
    color : tuple
        plot color
    name : pyqtgraph.LegendItem

    Returns
    -------
    plot : pyqtgraph.PlotItem
        2 x n plot with scatter plot and optional trace response plot for each feature (n)

    """
    if type(feature_list) is tuple:
        n_features = len(feature_list)
    else:
        n_features = 1
    if plot is None:
        plot = PlotGrid()
        plot.set_shape(n_features, 2)
        plot.show()
        for g in range(n_features):
            plot[g, 1].addLegend()

    for feature in range(n_features):
        if n_features > 1:
            current_feature = feature_list[feature]
            if median is True:
                mean = np.nanmedian(current_feature)
            else:
                mean = np.nanmean(current_feature)
            label = labels[feature]
            title = titles[feature]
        else:
            current_feature = feature_list
            mean = np.nanmean(current_feature)
            label = labels
            title = titles
        plot[feature, 0].setLabels(left=(label[0], label[1]))
        plot[feature, 0].hideAxis('bottom')
        plot[feature, 0].setTitle(title)
        if grand_trace is not None:
            plot[feature, 1].plot(grand_trace.time_values, grand_trace.data, pen=color, name=name)
        if len(current_feature) > 1:
            dx = pg.pseudoScatter(np.array(current_feature).astype(float), 0.7, bidir=True)
            plot[feature, 0].plot([i], [mean], symbol='o', symbolSize=20, symbolPen='k', symbolBrush=color)
            sem = stats.sem(current_feature, nan_policy='omit')
            if len(color) != 3:
                new_color = pg.glColor(color)
                color = (new_color[0]*255, new_color[1]*255, new_color[2]*255)
            plot[feature, 0].plot((0.3 * dx / dx.max()) + i, current_feature, pen=None, symbol='o', symbolSize=10, symbolPen='w',
                                symbolBrush=(color[0], color[1], color[2], 100))
        else:
            plot[feature, 0].plot([i], current_feature, pen=None, symbol='o', symbolSize=10, symbolPen='w',
                                symbolBrush=color)

    return plot


def get_expts(all_expts, cre_type, calcium=None, age=None, start=None, stop=None, dist_plot=None, color=None):
    """
    Collect and filter experiments by cre-type, caclium concentration, age, and date-range

    Parameters
    ----------
    all_expts : ExperimentList object
        list of experiments
    cre_type : string tuple
        cre-types of pre- and postsynaptic type of interest
    calcium : string
        can be 'high' or 'low' to filter for 2.0 mM external Ca++ or 1.3 mM
    age : string
        desired age range separated with '-'
    start : datetime
        start date of desired experiments
    stop : datetime
        stop date of desired experiments
    dist_plot : pyqtgraph.PlotItem
        plot to display distance profile of connections
    color : color tuple
        color for plot

    Returns
    -------
    expts : ExperimentList object
        list of filtered experiments
    legend : string
        string detailing cre-type, age, and calcium
    dist_plot : pyqtgraph.PlotItem
        plot of distance profile of connections from expts

    """
    expts = all_expts.select(cre_type=cre_type, calcium=calcium, age=age, start=start, stop=stop)
    if calcium is not None:
        if calcium == 'high':
            mM = '2.0 mM'
        elif calcium == 'low':
            mM = '1.3 mM'
    else:
        mM = ''
    if age is not None:
        age2 = 'P%s' % age
    else:
        age2 = 'All'
    legend = ("%s->%s, calcium = %s, age = %s" % (cre_type[0], cre_type[1], mM, age2))
    dist_plot = expts.distance_plot(cre_type[0], cre_type[1], plots=dist_plot, color=color, name=legend)[0]
    return expts, legend, dist_plot
