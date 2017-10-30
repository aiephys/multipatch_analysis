import pyqtgraph as pg
import numpy as np
import os
import pickle
from experiment_list import ExperimentList
from manuscript_figures import cache_response, get_amplitude, response_filter, trace_plot, bsub, write_cache, \
    induction_summary, recovery_summary, train_amp
from synapse_comparison import load_cache
from graphics import MatrixItem
from rep_connections import connections
from neuroanalysis.data import TraceList
from constants import INHIBITORY_CRE_TYPES, EXCITATORY_CRE_TYPES
from scipy import stats


pg.dbg()
app = pg.mkQApp()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

all_expts = ExperimentList(cache='expts_cache.pkl')
cre_types = ['sim1', 'tlx3', 'pvalb', 'sst', 'vip']
shape = (len(cre_types), len(cre_types))
calcium = 'high'
age = '40-60'
holding_e = [-68, -72]
holding_i = [-53, -60]
freqs = [10, 20, 50, 100, 200]
t_rec = [250, 500, 1000, 2000, 4000]
sweep_threshold = 5


pulse_cache_file = 'pulse_response_cache.pkl'
pulse_response_cache = load_cache(pulse_cache_file)
pulse_cache_change = []
train_cache_file = 'train_response_cache.pkl'
train_response_cache = load_cache(train_cache_file)
train_cache_change = []

e_plot = pg.GraphicsWindow()
e_row = 0
i_plot = pg.GraphicsWindow()
i_row = 0

for c1, pre_type in enumerate(cre_types):
    for c2, post_type in enumerate(cre_types):
        grand_pulse_response = []
        grand_induction = {}
        grand_recovery = {}
        offset_ind = {}
        offset_rec = {}
        expt_list = all_expts.select(cre_type=[pre_type, post_type], calcium=calcium, age=age)
        if pre_type in EXCITATORY_CRE_TYPES:
            p1 = e_plot.addPlot(row=e_row, col=0)
            p2 = e_plot.addPlot(row=e_row, col=1)
            p3 = e_plot.addPlot(row=e_row, col=2)
            p4 = e_plot.addPlot(row=e_row, col=3)
            p5 = e_plot.addPlot(row=e_row, col=4)
            holding = holding_e
            sign = '+'
            avg_color = {'color': 'r', 'width': 2}
            e_row += 1
        else:
            p1 = i_plot.addPlot(row=i_row, col=0)
            p2 = i_plot.addPlot(row=i_row, col=1)
            p3 = i_plot.addPlot(row=i_row, col=2)
            p4 = i_plot.addPlot(row=i_row, col=3)
            p5 = i_plot.addPlot(row=i_row, col=4)
            holding = holding_i
            sign = '-'
            avg_color = {'color': 'b', 'width': 2}
            i_row += 1
        if e_row == 0 or i_row:
            p1.setTitle('First Pulse Response, Example Connection')
            p2.setTitle('First Pulse Response, All Connections')
            p3.setTitle('50 Hz Train Response')
            p4.setTitle('Induction Amplitudes')
            p5.setTitle('Recovery Amplitudes')
        trace_color = (0, 0, 0, 30)
        trace_color2 = (255, 0, 255, 30)
        p1.addLegend()
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
                        avg_trace, _, amp_sign, _ = get_amplitude(pulse_subset)
                        if pre_type in EXCITATORY_CRE_TYPES and amp_sign is '-':
                            continue
                        elif pre_type in INHIBITORY_CRE_TYPES and amp_sign is '+':
                            continue
                        avg_trace.t0 = 0
                        grand_pulse_response.append(avg_trace)
                        if [expt.uid, pre, post] == connections[pre_type, post_type]:
                            for sweep in pulse_subset:
                                bsub_sweep = bsub(sweep)
                                p1 = trace_plot(bsub_sweep, color=trace_color, plot=p1, x_range=[0, 27e-3])
                            p1 = trace_plot(avg_trace, color=(255, 0, 255), plot=p1, x_range=[0, 27e-3],
                                            name=('%s -> %s' % (pre_type, post_type)))
                            p2 = trace_plot(avg_trace, color=trace_color2, plot=p2, x_range=[0, 27e-3])
                        else:
                           p2 = trace_plot(avg_trace, color=trace_color, plot=p2, x_range=[0, 27e-3])

                    train_response, cache_change = cache_response(expt, pre, post, train_response_cache, type='train')
                    train_cache_change.append(cache_change)
                    grand_induction, offset_ind = induction_summary(train_response, freqs, holding, thresh=sweep_threshold,
                                                        ind_dict=grand_induction, offset_dict=offset_ind)
                    grand_recovery, offset_rec = recovery_summary(train_response, t_rec, holding, thresh=sweep_threshold,
                                                      rec_dict=grand_recovery, offset_dict=offset_rec)

        if len(grand_pulse_response) > 0:
            grand_pulse_trace = TraceList(grand_pulse_response).mean()
            p2 = trace_plot(grand_pulse_trace, color=avg_color, plot=p2, x_range=[0, 27e-3])
            if len(grand_induction) > 0:
                for f, freq in freqs:
                    if freq == 50:
                        grand_ind_trace = TraceList(grand_induction[freq][0]).mean()
                        grand_rec_trace = TraceList(grand_induction[freq][1]).mean()
                        for ind in grand_induction[freq][0]:
                            p3 = trace_plot(ind, color=trace_color, plot=p3)
                        for rec in grand_induction[freq][1]:
                            p3 = trace_plot(rec, color=trace_color, plot=p3)
                        p3 = trace_plot(grand_ind_trace, color=avg_color, plot=p3)
                        p3 = trace_plot(grand_rec_trace, color=avg_color, plot=p3)
                    offset = offset_ind[freq]
                    ind_amp = train_amp(grand_induction[freq], offset, sign)
                    grand_ind_amp = np.mean(ind_amp, 0)
                    ind_amp_sem = stats.sem(ind_amp)
                    f_color = pg.hsvColor(1, sat=float(f+0.5)/len(freqs), val=1)
                    p4.plot(grand_ind_amp/grand_ind_amp[0], name=('  %d Hz' % freq), pen=f_color, symbol='t',
                            symbolBrush=f_color, symbolPen=None)
            if len(grand_recovery) > 0:
                for t, delta in t_rec:
                    offset = offset_rec[delta]
                    rec_amp = train_amp(grand_recovery, offset, sign)
                    grand_rec_amp = np.mean(rec_amp)
                    grand_rec_sem = stats.sem(rec_amp)
                    t_color = pg.hsvColor(1, sat=float(t+0.5)/len(t_rec), val=1)
                    p5.plot(grand_rec_amp/grand_rec_amp[0], name=('  %d ms' % delta), pen=t_color, symbol='t',
                            symbolBrush=t_color, symbolPen=None)

    if sum(pulse_cache_change) > 0:
        write_cache(pulse_response_cache, pulse_cache_file)
    if sum(train_cache_change) > 0:
        write_cache(train_response_cache, train_cache_file)

