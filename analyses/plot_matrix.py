import pyqtgraph as pg
import numpy as np
import colorsys as colorsys
import os
import pickle
from experiment_list import ExperimentList
from manuscript_figures import cache_response, get_amplitude, response_filter, trace_plot, bsub, write_cache, \
    induction_summary, recovery_summary, train_amp, pulse_qc, train_qc, subplots
from synapse_comparison import load_cache
from aisynphys.ui.graphics import MatrixItem
from rep_connections import connections
from neuroanalysis.data import TSeriesList
from aisynphys.constants import INHIBITORY_CRE_TYPES, EXCITATORY_CRE_TYPES
from aisynphys.experiment_list import cached_experiments
from scipy import stats


pg.dbg()
app = pg.mkQApp()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

all_expts = cached_experiments()
cre_types = ['sim1', 'tlx3', 'pvalb', 'sst', 'vip']
no_connections = [('sim1', 'tlx3'), ('tlx3', 'sim1'), ('sim1', 'vip'), ('tlx3', 'vip'), ('sst', 'sim1'), ('vip', 'sim1'),
                  ('vip', 'tlx3'), ('vip', 'pvalb'), ('sim1', 'sst')]
shape = (len(cre_types), len(cre_types))
calcium = 'high'
age = '40-60'
holding_e = [-68, -72]
holding_i = [-53, -60]
freqs = [10, 20, 50, 100, 200]
t_rec = [250, 500, 1000, 2000, 4000]
sweep_threshold = 5
amp_thresh = 100e-6
amp = None

plt = pg.plot()

pulse_cache_file = 'pulse_response_cache.pkl'
pulse_response_cache = load_cache(pulse_cache_file)
pulse_cache_change = []
#train_cache_file = 'train_response_cache.pkl'
#train_response_cache = load_cache(train_cache_file)
train_cache_change = []

big_plot = pg.GraphicsLayoutWidget()
big_plot.show()
row = 0

pulse_amp = {}
ind_index = {}
rec_index = {}

for c1, pre_type in enumerate(cre_types):
    for c2, post_type in enumerate(cre_types):
        if (pre_type, post_type) in no_connections:
            continue
        p1, p2, p3, p4, p5 = subplots(name=big_plot, row=row)
        key = (pre_type, post_type)
        train_cache_file = ('%s-%s_train_response.pkl' % (pre_type, post_type))
        train_response_cache = load_cache(train_cache_file)
        grand_pulse_response = []
        grand_induction = {}
        grand_recovery = {}
        offset_ind = {}
        offset_rec = {}
        expt_list = all_expts.select(cre_type=[pre_type, post_type], calcium=calcium, age=age)
        if pre_type in EXCITATORY_CRE_TYPES and post_type in EXCITATORY_CRE_TYPES:
            holding = holding_e
            sign = '+'
            avg_color = {'color': (255, 0, 0), 'width': 2}
            hue = colorsys.rgb_to_hsv(255, 0, 0)
            hue = hue[0]
        elif pre_type in EXCITATORY_CRE_TYPES and post_type in INHIBITORY_CRE_TYPES:
            holding = holding_i
            sign = '+'
            avg_color = {'color': (255, 140, 0), 'width': 2}
            hue = colorsys.rgb_to_hsv(255, 140, 0)
            hue = hue[0]
        elif pre_type in INHIBITORY_CRE_TYPES and post_type in EXCITATORY_CRE_TYPES:
            holding = holding_i
            sign = '-'
            avg_color = {'color': (138, 43, 226), 'width': 2}
            hue = colorsys.rgb_to_hsv(138, 43, 226)
            hue = hue[0]
        elif pre_type in INHIBITORY_CRE_TYPES and post_type in INHIBITORY_CRE_TYPES:
            holding = holding_i
            sign = '-'
            avg_color = {'color': (0, 0, 255), 'width': 2}
            hue = colorsys.rgb_to_hsv(0, 0, 255)
            hue = hue[0]
        trace_color = (0, 0, 0, 30)
        trace_color2 = (255, 0, 255, 30)
        p1.addLegend()
        p2.addLegend()
        p3.addLegend()
        p4.addLegend()
        p4.setLabels(left=('Norm Amp', ''), bottom=('Pulse Number', ''))
        p5.addLegend()
        p5.setLabels(left=('Norm Amp', ''), bottom=('Pulse Number', ''))
        for expt in expt_list:
            for pre, post in expt.connections:
                if expt.cells[pre].cre_type == pre_type and expt.cells[post].cre_type == post_type:
                    pulse_response, cache_change = cache_response(expt, pre, post, pulse_response_cache, type='pulse')
                    pulse_cache_change.append(cache_change)
                    pulse_subset = response_filter(pulse_response, freq_range=[0, 50], holding_range=holding, pulse=True)
                    if len(pulse_subset) >= sweep_threshold:
                        plt.clear()
                        pass_qc = pulse_qc(pulse_subset, baseline=4, pulse=4, plot=plt)
                        if len(pass_qc) >= sweep_threshold:
                            avg_trace, amp, amp_sign, _ = get_amplitude(pass_qc)
                            if pre_type in EXCITATORY_CRE_TYPES and amp_sign is '-':
                                continue
                            elif pre_type in INHIBITORY_CRE_TYPES and amp_sign is '+':
                                continue
                            if key not in pulse_amp.keys():
                                pulse_amp[key] = []
                            pulse_amp[key].append(amp)
                            avg_trace.t0 = 0
                            grand_pulse_response.append(avg_trace)
                            if [expt.uid, pre, post] == connections[pre_type, post_type]:
                                for sweep in pass_qc:
                                    p1 = trace_plot(sweep, color=trace_color, plot=p1, x_range=[0, 27e-3])
                                p1 = trace_plot(avg_trace, color=(255, 0, 255), plot=p1, x_range=[0, 27e-3],
                                                name=('%s -> %s' % (pre_type, post_type)))
                                p2 = trace_plot(avg_trace, color=trace_color2, plot=p2, x_range=[0, 27e-3])
                            else:
                               p2 = trace_plot(avg_trace, color=trace_color, plot=p2, x_range=[0, 27e-3])

                    if amp is not None and abs(amp) > amp_thresh:
                        train_response, cache_change = cache_response(expt, pre, post, train_response_cache, type='train')
                        train_cache_change.append(cache_change)
                        if (50, 0.25) in [(k[0], np.round(k[1], 2)) for k in train_response['responses'].keys()]:
                            grand_induction, offset_ind = induction_summary(train_response, freqs, holding, thresh=sweep_threshold,
                                                                ind_dict=grand_induction, offset_dict=offset_ind)
                            grand_recovery, offset_rec = recovery_summary(train_response, t_rec, holding, thresh=sweep_threshold,
                                                        rec_dict=grand_recovery, offset_dict=offset_rec)

        if len(grand_pulse_response) > 0:
            grand_pulse_trace = TSeriesList(grand_pulse_response).mean()
            p2 = trace_plot(grand_pulse_trace, color=avg_color, plot=p2, x_range=[0, 27e-3], name=('n = %d' % len(grand_pulse_response)))
            if len(grand_induction) > 0:
                for f, freq in enumerate(freqs):
                    if freq in grand_induction:
                        offset = offset_ind[freq]
                        ind_pass_qc = train_qc(grand_induction[freq], offset, amp=amp_thresh, sign=sign)
                        n = len(ind_pass_qc[0])
                        if n > 0:
                            ind_amp = train_amp(ind_pass_qc, offset, sign)
                            grand_ind_amp = np.nanmean(ind_amp, 0)
                            ind_amp_sem = stats.sem(ind_amp)
                            if freq not in ind_index.keys():
                                ind_index[freq] = {}
                            if key not in ind_index[freq].keys():
                                ind_index[freq][key] = []
                            for n in range(ind_amp.shape[0]):
                                ind_index[freq][key].append(ind_amp[n, 7] / ind_amp[n, 0])
                            if freq == 50:
                                grand_ind_trace = TSeriesList(ind_pass_qc[0]).mean()
                                grand_rec_trace = TSeriesList(ind_pass_qc[1]).mean()
                                for ind in ind_pass_qc[0]:
                                    p3 = trace_plot(ind, color=trace_color, plot=p3)
                                for rec in ind_pass_qc[1]:
                                    p3 = trace_plot(rec, color=trace_color, plot=p3)
                                p3 = trace_plot(grand_ind_trace, color=avg_color, plot=p3)
                                p3 = trace_plot(grand_rec_trace, color=avg_color, plot=p3, name=('n = %d'% n))

                            f_color = pg.hsvColor(hue=hue, sat=float(f+0.5)/len(freqs), val=1)
                            p4.plot(grand_ind_amp/grand_ind_amp[0], name=('  %d Hz, n = %d' % (freq, n)), pen=f_color, symbol='t',
                                    symbolBrush=f_color, symbolPen=None)
                            # ind_err = pg.ErrorBarItem(x=np.arange(12), y=np.array(grand_ind_amp/grand_ind_amp[0]),
                            #                       height=np.array(ind_amp_sem), beam=0.3)
                            # p4.addItem(ind_err)
            if len(grand_recovery) > 0:
                for t, delta in enumerate(t_rec):
                    if delta in grand_recovery:
                        offset = offset_rec[delta]
                        rec_pass_qc = train_qc(grand_recovery[delta], offset, amp=amp_thresh, sign=sign)
                        n = len(rec_pass_qc[1])
                        if n > 0:
                            rec_amp = train_amp(rec_pass_qc, offset, sign)
                            grand_rec_amp = np.mean(rec_amp, 0)
                            rec_amp_sem = stats.sem(rec_amp)
                            if delta not in rec_index.keys():
                                rec_index[delta] = {}
                            if key not in rec_index[delta].keys():
                                rec_index[delta][key] = []
                            for n in range(rec_amp.shape[0]):
                                rec_index[delta][key].append(rec_amp[n, 8] / rec_amp[n, 0])
                            t_color = pg.hsvColor(hue, sat=float(t+0.5)/len(t_rec), val=1)
                            p5.plot(grand_rec_amp/grand_rec_amp[0], name=('  %d ms, n = %d' % (delta, n)), pen=t_color, symbol='t',
                                    symbolBrush=t_color, symbolPen=None)
                            # rec_err = pg.ErrorBarItem(x=np.arange(12),y=np.array(grand_rec_amp / grand_rec_amp[0]),
                            #                           height=np.array(rec_amp_sem), beam=0.3)
                            # p5.addItem(rec_err)
        row += 1
    # if sum(pulse_cache_change) > 0:
    #     write_cache(pulse_response_cache, pulse_cache_file)
    # if sum(train_cache_change) > 0:
    #     write_cache(train_response_cache, train_cache_file)

feature_cache = {}
feature_cache['Amplitudes'] = pulse_amp
feature_cache['Induction'] = ind_index
feature_cache['Recovery'] = rec_index
write_cache(feature_cache, 'feature_cache_file.pkl')

