import pyqtgraph as pg
import numpy as np
import os
import pickle
from experiment_list import ExperimentList
from manuscript_figures import cache_response, get_amplitude, response_filter, trace_plot, bsub, trace_avg, induction_summary, recovery_summary
from synapse_comparison import load_cache
from graphics import MatrixItem
from rep_connections import connections
from neuroanalysis.data import TraceList
from constants import INHIBITORY_CRE_TYPES, EXCITATORY_CRE_TYPES


#pg.dbg()
app = pg.mkQApp()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

all_expts = ExperimentList(cache='expts_cache.pkl')
cre_types = ['sim1', 'tlx3', 'pvalb', 'sst', 'vip']
shape = (len(cre_types), len(cre_types))
calcium = 'high'
age = '40-51'
holding_e = [-68, -72]
holding_i = [-53, -60]
freqs = [10, 20, 50, 100, 200]
t_rec = [250, 500, 1000, 2000, 4000]
sweep_threshold = 5


pulse_cache_file = 'pulse_response_cache.pkl'
pulse_response_cache = load_cache(pulse_cache_file)
train_cache_file = 'D:\\train_response_cache.pkl'
train_response_cache = load_cache(train_cache_file)

e_plot = pg.GraphicsLayout()
i_plot = pg.GraphicsLayout()

for c1, pre_type in enumerate(cre_types):
    for c2, post_type in enumerate(cre_types):
        grand_pulse_response = []
        grand_induction = {}
        grand_recovery = {}
        expt_list = all_expts.select(cre_type=[pre_type, post_type], calcium=calcium, age=age)
        if c1 in EXCITATORY_CRE_TYPES:
            p1 = e_plot.addPlot(row=e_plot.nextRow(), col=0)
            p2 = e_plot.addPlot(row=e_plot.nextRow(), col=1)
            p3 = e_plot.addPlot(row=e_plot.nextRow(), col=2)
            #p4 = e_plot.addPlot(row=e_plot.nextRow(), col=3)
        else:
            p1 = e_plot.addPlot(row=i_plot.nextRow(), col=0)
            p2 = e_plot.addPlot(row=i_plot.nextRow(), col=1)
            p3 = e_plot.addPlot(row=i_plot.nextRow(), col=2)
            # p4 = e_plot.addPlot(row=i_plot.nextRow(), col=3)
        trace_color = (0, 0, 0, 30)
        trace_color2 = (255, 0, 255, 30)
        avg_color = {'color': (50, 205, 50), 'width': 2}
        for expt in expt_list:
            for pre, post in expt.connections:
                if expt.cells[pre].cre_type == pre_type and expt.cells[post].cre_type == post_type:
                    pulse_response = cache_response(expt, pre, post, pulse_cache_file, pulse_response_cache, type='pulse')
                    if pre_type in EXCITATORY_CRE_TYPES:
                       holding = holding_e
                    elif pre_type in INHIBITORY_CRE_TYPES:
                       holding = holding_i
                    pulse_subset = response_filter(pulse_response, freq_range=[0, 50], holding_range=holding)
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
                            p1 = trace_plot(avg_trace, color=(255, 0, 255), plot=p1, x_range=[0, 27e-3])
                            p2 = trace_plot(avg_trace, color=trace_color2, plot=p2, x_range=[0, 27e-3])
                        else:
                           p2 = trace_plot(avg_trace, color=trace_color, plot=p2, x_range=[0, 27e-3])

                    train_response = cache_response(expt, pre, post,train_cache_file, train_response_cache, type='train')
                    grand_induction = induction_summary(train_response, freqs, holding, thresh=sweep_threshold, ind_dict=grand_induction)
                    grand_recovery = recovery_summary(train_response, t_rec, holding, thresh=sweep_threshold, rec_dict=grand_recovery)

        if len(grand_pulse_response) > 0:
            grand_pulse_trace = TraceList(grand_pulse_response).mean()
            p2 = trace_plot(grand_pulse_trace, color=avg_color, plot=p2, x_range=[0, 27e-3])
            for f in freqs:
                if f == 50:
                    grand_ind_trace = TraceList(grand_induction[f][0]).mean()
                    grand_rec_trace = TraceList(grand_recovery[f][0]).mean()
                    p3 = [trace_plot(ind, color=trace_color, plot=p3) for ind in grand_induction[f][0]]
                    p3 = trace_plot(grand_ind_trace, color=avg_color, plot=p3)
                    p3 = trace_plot(grand_rec_trace, color=avg_color, plot=p3)

    pickle.dump(pulse_response_cache, open(pulse_cache_file+'.new', 'wb'))
    if os.path.exists(pulse_cache_file):
        os.path.remove(pulse_cache_file)
    os.rename(pulse_cache_file+'.new', pulse_cache_file)

    pickle.dump(train_response_cache, open(train_cache_file + '.new', 'wb'))
    if os.path.exists(train_cache_file):
        os.path.remove(train_cache_file)
    os.rename(train_cache_file + '.new', train_cache_file)

