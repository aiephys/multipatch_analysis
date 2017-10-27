import pyqtgraph as pg
import numpy as np
import os
import pickle
from experiment_list import ExperimentList
from manuscript_figures import cache_response, get_amplitude, response_filter, trace_plot, bsub, trace_avg
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
sweep_threshold = 5


pulse_cache_file = 'pulse_response_cache.pkl'
pulse_response_cache = load_cache(pulse_cache_file)
train_cache_file = 'D:\\train_response_cache.pkl'
train_response_cache = load_cache(train_cache_file)

win = pg.GraphicsLayoutWidget()
win.show()
vb = win.addViewBox()
vb.invertY()
vb.setAspectLocked()
expts = all_expts.select(calcium=calcium, age=age)
bgcolor = expts.matrix(cre_types, cre_types)
fgcolor = np.empty(shape, dtype='U')
fgcolor[:] = 'w'
text = np.empty(shape, dtype='U')
text[:] = ''
matrix = MatrixItem(text=text, fgcolor=fgcolor, bgcolor=bgcolor, rows=cre_types, cols=cre_types, size=400)
vb.addItem(matrix)
element_bgcolor = (255, 255, 255, 200)

for c1, pre_type in enumerate(cre_types):
    row_label = matrix.row_labels[c1]
    row_label.setDefaultTextColor(pg.mkColor('k'))
    row_label.setTextWidth(400)
    for c2, post_type in enumerate(cre_types):
        grand_pulse_response = []
        grand_train_response = {'ind': [], 'rec': []}
        expt_list = all_expts.select(cre_type=[pre_type, post_type], calcium=calcium, age=age)
        col_label = matrix.col_labels[c2]
        col_label.setDefaultTextColor(pg.mkColor('k'))
        col_label.setTextWidth(400)
        element = pg.GraphicsLayout()
        p1 = element.addPlot(row=0, col=0, rowspan=1, colspan=1)
        p2 = element.addPlot(row=0, col=1, rowspan=1, colspan=1)
        p3 = element.addPlot(row=1, col=0, rowspan=1, colspan=2)
        #p4 = element.addPlot(row=1, col=1, rowspan=1, colspan=1)
        element.setParentItem(matrix.cells[c1][c2])
        element.resize(400, 400)
        trace_color = (0, 0, 0, 30)
        trace_color2 = (255, 0, 255, 30)
        avg_color = {'color': (50, 205, 50), 'width': 2}
        for expt in expt_list:
            for pre, post in expt.connections:
                if expt.cells[pre].cre_type == pre_type and expt.cells[post].cre_type == post_type:
                    pulse_response = cache_response(expt, pre, post, pulse_response_cache, type='pulse')
                    # if pre_type in EXCITATORY_CRE_TYPES:
                    #    holding = holding_e
                    # elif pre_type in INHIBITORY_CRE_TYPES:
                    #    holding = holding_i
                    # pulse_subset = response_filter(pulse_response, freq_range=[0, 50], holding_range=holding)
                    # if len(pulse_subset) >= sweep_threshold:
                    #     avg_trace, _, amp_sign, _ = get_amplitude(pulse_subset)
                    #     if pre_type in EXCITATORY_CRE_TYPES and amp_sign is '-':
                    #         continue
                    #     elif pre_type in INHIBITORY_CRE_TYPES and amp_sign is '+':
                    #         continue
                    #     avg_trace.t0 = 0
                    #     grand_pulse_response.append(avg_trace)
                    #     p1.vb.setBackgroundColor(element_bgcolor)
                    #     p2.vb.setBackgroundColor(element_bgcolor)
                    #     if [expt.uid, pre, post] == connections[pre_type, post_type]:
                    #         for sweep in pulse_subset:
                    #             bsub_sweep = bsub(sweep)
                    #             p1 = trace_plot(bsub_sweep, color=trace_color, plot=p1, x_range=[0, 27e-3])
                    #         p1 = trace_plot(avg_trace, color=(255, 0, 255), plot=p1, x_range=[0, 27e-3])
                    #         p2 = trace_plot(avg_trace, color=trace_color2, plot=p2, x_range=[0, 27e-3])
                    #     else:
                    #        p2 = trace_plot(avg_trace, color=trace_color, plot=p2, x_range=[0, 27e-3])

                    train_response = cache_response(expt, pre, post, train_response_cache, type='train')
        #             induction = response_filter(train_response[0], freq_range=[50, 50], holding_range=holding)
        #             recovery = response_filter(train_response[1], freq_range=[50, 50], holding_range=holding, delta_t=250)
        #             if len(recovery) >= sweep_threshold:
        #                 ind_avg = trace_avg(induction)
        #                 rec_avg = trace_avg(recovery)
        #                 rec_avg.t0 = ind_avg.time_values[-1] + 0.1
        #                 grand_train_response['ind'].append(ind_avg)
        #                 grand_train_response['rec'].append(rec_avg)
        #                 p3.vb.setBackgroundColor(element_bgcolor)
        #                 #p4.vb.setBackgroundColor(element_bgcolor)
        #                 if [expt.uid, pre, post] == connections[pre_type, post_type]:
        #                     p3 = trace_plot(ind_avg, color=trace_color2, plot=p3)
        #                     p3 = trace_plot(rec_avg, color=trace_color2, plot=p3)
        #                 else:
        #                     p3 = trace_plot(ind_avg, color=trace_color, plot=p3)
        #                     p3 = trace_plot(rec_avg, color=trace_color, plot=p3)
        #         # else:
        #         #     element.removeItem(p1)
        #         #     element.removeItem(p2)
        #         #     element.removeItem(p3)
        #         #     element.removeItem(p4)
        #
        # if len(grand_pulse_response) > 0:
        #     grand_pulse_trace = TraceList(grand_pulse_response).mean()
        #     p2 = trace_plot(grand_pulse_trace, color=avg_color, plot=p2, x_range=[0, 27e-3])
        #     grand_ind_trace = TraceList(grand_train_response['ind']).mean()
        #     grand_rec_trace = TraceList(grand_train_response['rec']).mean()
        #     p3 = trace_plot(grand_ind_trace, color=avg_color, plot=p3)
        #     p3 = trace_plot(grand_rec_trace, color=avg_color, plot=p3)

    pickle.dump(pulse_response_cache, open(pulse_cache_file+'.new', 'wb'))
    if os.path.exists(pulse_cache_file):
        os.path.remove(pulse_cache_file)
    os.rename(pulse_cache_file+'.new', pulse_cache_file)

    pickle.dump(train_response_cache, open(train_cache_file + '.new', 'wb'))
    if os.path.exists(train_cache_file):
        os.path.remove(train_cache_file)
    os.rename(train_cache_file + '.new', train_cache_file)