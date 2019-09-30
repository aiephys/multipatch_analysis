import pyqtgraph as pg
import sys
import argparse
import numpy as np
import colorsys
from aisynphys.experiment_list import cached_experiments
from manuscript_figures import cache_response, train_amp, induction_summary, recovery_summary, write_cache, get_response, \
    train_qc, colors_human, colors_mouse, deconv_train
from synapse_comparison import load_cache
from neuroanalysis.data import TSeriesList
from neuroanalysis.ui.plot_grid import PlotGrid
from aisynphys.synaptic_dynamics import RawDynamicsAnalyzer
from rep_connections import ee_connections, human_connections, all_connections, ii_connections, ei_connections, ie_connections, no_include
pg.dbg()
app = pg.mkQApp()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

parser = argparse.ArgumentParser(description='Enter organism and type of connection you"d like to analyze ex: mouse ee (all mouse excitatory-'
                'excitatory). Alternatively enter a cre-type connection ex: sim1-sim1')
parser.add_argument('--organism', dest='organism', help='Select mouse or human')
parser.add_argument('--connection', dest='connection', help='Specify connections to analyze')
args = vars(parser.parse_args(sys.argv[1:]))

all_expts = cached_experiments()

if args['organism'] == 'mouse':
    color_palette = colors_mouse
    calcium = 'high'
    age = '40-60'
    connection = args['connection']
    if connection == 'ee':
        connection_types = ee_connections.keys()
        sign = '+'
        threshold = [None, None, 0.03e-3, 0.03e-3, None]
        qc_params = (sign, [0.5e-3, 1e-3, 1e-3, 1e-3, 1e-3])
    elif connection =='ii':
        connection_types = ii_connections.keys()
        sign = '-'
    elif connection == 'ei':
        connection_types = ei_connections.keys()
        sign = '+'
    elif connection == 'ie':
        connection_types = ie_connections.keys()
        sign = '-'
    elif connection == 'all':
        connection_types = all_connections.keys()
        sign = None
    elif len(connection.split('-')) == 2:
        c_type = connection.split('-')
        if c_type[0] == '2/3':
            pre_type = ('2/3', 'unknown')
        else:
            pre_type = (None, c_type[0])
        if c_type[1] == '2/3':
            post_type = ('2/3', 'unknown')
        else:
            post_type = (None, c_type[0])
        connection_types = [(pre_type, post_type)]
        threshold = [None]
        qc_params = ('+', [1e-3])
elif args['organism'] == 'human':
    color_palette = colors_human
    calcium = None
    age = None
    connection = args['connection']
    if connection == 'ee':
        connection_types = human_connections.keys()
        sign = '+'
        threshold = [None, None, None, None]
        qc_params = (sign, [1e-3, 1e-3, 1e-3, 1e-3])
    else:
        c_type = connection.split('-')
        connection_types = [((c_type[0], 'unknown'), (c_type[1], 'unknown'))]
        sign='+'
        threshold = [None]
        qc_params = (sign, [1e-3])

holding = [-55, -75]
freqs = [10, 20, 50, 100]
rec_t = [250, 500, 1000, 2000, 4000]
sweep_threshold = 3
deconv = True

# cache_file = 'train_response_cache.pkl'
# response_cache = load_cache(cache_file)
# cache_change = []
# log_rec_plt = pg.plot()
# log_rec_plt.setLogMode(x=True)
qc_plot = pg.plot()
ind_plot = PlotGrid()
ind_plot.set_shape(4, len(connection_types))
ind_plot.show()
rec_plot = PlotGrid()
rec_plot.set_shape(5, len(connection_types))
rec_plot.show()
if deconv is True:
    deconv_ind_plot = PlotGrid()
    deconv_ind_plot.set_shape(4, len(connection_types))
    deconv_ind_plot.show()
    deconv_rec_plot = PlotGrid()
    deconv_rec_plot.set_shape(5,len(connection_types))
    deconv_rec_plot.show()
summary_plot = PlotGrid()
summary_plot.set_shape(len(connection_types), 2)
summary_plot.show()
symbols = ['o', 's', 'd', '+', 't']
trace_color = (0, 0, 0, 5)
ind_amp_summary = {}
ind_uid = {}
rec_amp_summary = {}
rec_uid = {}
for c in range(len(connection_types)):
    induction_grand = {}
    recovery_grand = {}
    pulse_offset_ind = {}
    pulse_offset_rec = {}
    cre_type = (connection_types[c][0][1], connection_types[c][1][1])
    target_layer = (connection_types[c][0][0], connection_types[c][1][0])
    c_type = connection_types[c]
    expt_list = all_expts.select(cre_type=cre_type, target_layer=target_layer, calcium=calcium, age=age)
    color = color_palette[c]
    hue = colorsys.rgb_to_hsv(color[0], color[1], color[2])
    hue = hue[0]
    summary_plot[c, 0].addLegend()
    summary_plot[c, 1].addLegend()
    for expt in expt_list:
        if expt.connections is None:
            continue
        for pre, post in expt.connections:
            if [expt.uid, pre, post] in no_include:
                continue
            if expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]:
                train_response, artifact = get_response(expt, pre, post, type='train')
                if threshold[c] is not None and artifact > threshold[c]:
                    continue
                induction_grand, pulse_offset_ind = induction_summary(train_response, freqs, holding, thresh=sweep_threshold,
                                                                ind_dict=induction_grand, offset_dict=pulse_offset_ind,
                                                                      uid=(expt.uid, pre, post))
                recovery_grand, pulse_offset_rec = recovery_summary(train_response, rec_t, holding, thresh=sweep_threshold,
                                                                rec_dict=recovery_grand, offset_dict=pulse_offset_rec,
                                                                    uid=(expt.uid, pre, post))
    for f, freq in enumerate(freqs):
        if freq not in induction_grand.keys():
            print ("%d Hz not represented in data set for %s" % (freq, c_type))
            continue
        ind_offsets = pulse_offset_ind[freq]
        qc_plot.clear()
        ind_pass_qc = train_qc(induction_grand[freq], ind_offsets, amp=qc_params[1][c], sign=qc_params[0], plot=qc_plot)
        n_synapses = len(ind_pass_qc[0])
        if n_synapses > 0:
            induction_grand_trace = TSeriesList(ind_pass_qc[0]).mean()
            ind_rec_grand_trace = TSeriesList(ind_pass_qc[1]).mean()
            ind_amp = train_amp(ind_pass_qc, ind_offsets, '+')
            ind_amp_grand = np.nanmean(ind_amp, 0)

            if f == 0:
                ind_plot[f, c].setTitle(connection_types[c])
                type = pg.LabelItem('%s -> %s' % connection_types[c])
                type.setParentItem(summary_plot[c, 0])
                type.setPos(50, 0)
            if c == 0:
                label = pg.LabelItem('%d Hz Induction' % freq)
                label.setParentItem(ind_plot[f, c].vb)
                label.setPos(50, 0)
                summary_plot[c, 0].setTitle('Induction')
            ind_plot[f, c].addLegend()
            [ind_plot[f, c].plot(ind.time_values, ind.data, pen=trace_color) for ind in ind_pass_qc[0]]
            [ind_plot[f, c].plot(rec.time_values, rec.data, pen=trace_color) for rec in ind_pass_qc[1]]

            ind_plot[f, c].plot(induction_grand_trace.time_values, induction_grand_trace.data, pen={'color': color, 'width': 2},
                                name=("n = %d" % n_synapses))
            ind_plot[f, c].plot(ind_rec_grand_trace.time_values, ind_rec_grand_trace.data, pen={'color': color, 'width': 2})
            ind_plot[f, c].setLabels(left=('Vm', 'V'))
            ind_plot[f, c].setLabels(bottom=('t', 's'))
            if deconv is True:
                ind_deconv = deconv_train(ind_pass_qc[:2])
                ind_deconv_grand = TSeriesList(ind_deconv[0]).mean()
                ind_rec_deconv_grand = TSeriesList(ind_deconv[1]).mean()
                [deconv_ind_plot[f, c].plot(ind.time_values, ind.data, pen=trace_color) for ind in ind_deconv[0]]
                [deconv_ind_plot[f, c].plot(rec.time_values, rec.data, pen=trace_color) for rec in ind_deconv[1]]
                deconv_ind_plot[f, c].plot(ind_deconv_grand.time_values, ind_deconv_grand.data,
                                    pen={'color': color, 'width': 2}, name=("n = %d" % n_synapses))
                deconv_ind_plot[f, c].plot(ind_rec_deconv_grand.time_values, ind_rec_deconv_grand.data,
                                pen={'color': color, 'width': 2})
            summary_plot[c, 0].setLabels(left=('Norm Amp', ''))
            summary_plot[c, 0].setLabels(bottom=('Pulse Number', ''))
            f_color = pg.hsvColor(hue=hue, sat=float(f+0.5)/len(freqs), val=1)
            summary_plot[c, 0].plot(ind_amp_grand/ind_amp_grand[0], name=('  %d Hz' % freq), pen=f_color, symbol=symbols[f],
                                symbolBrush=f_color, symbolPen=None)
            if connection_types[c] not in ind_amp_summary.keys():
                ind_amp_summary[connection_types[c]] = []
                ind_uid[connection_types[c]] = []
            ind_amp_summary[connection_types[c]].append([freq, ind_amp])
            ind_uid[connection_types[c]].append([freq, ind_pass_qc[2]])

    for t, delta in enumerate(rec_t):
        if delta not in recovery_grand.keys():
            print ("%d ms not represented in data set for %s" % (delta, c_type))
            continue
        rec_offsets = pulse_offset_rec[delta]
        rec_pass_qc = train_qc(recovery_grand[delta], rec_offsets, amp=qc_params[1][c], sign=qc_params[0], plot=None)
        n_synapses = len(rec_pass_qc[0])
        if n_synapses > 0:
            recovery_grand_trace = TSeriesList(rec_pass_qc[1]).mean()
            rec_ind_grand_trace = TSeriesList(rec_pass_qc[0]).mean()
            rec_amp = train_amp(rec_pass_qc, rec_offsets, '+')
            rec_amp_grand = np.mean(rec_amp, 0)
            if t == 0:
                rec_plot[t, c].setTitle(connection_types[c])
            if c == 0:
                label = pg.LabelItem('%d ms Recovery' % delta)
                label.setParentItem(rec_plot[t, c].vb)
                label.setPos(50, 0)
                summary_plot[c, 1].setTitle('Recovery')
            rec_plot[t, c].addLegend()
            [rec_plot[t, c].plot(ind.time_values, ind.data, pen=trace_color) for ind in rec_pass_qc[0]]
            [rec_plot[t, c].plot(rec.time_values, rec.data, pen=trace_color) for rec in rec_pass_qc[1]]
            rec_plot[t, c].plot(rec_ind_grand_trace.time_values, rec_ind_grand_trace.data, pen={'color': color, 'width': 2},
                                name=("n = %d" % n_synapses))
            rec_plot[t, c].plot(recovery_grand_trace.time_values, recovery_grand_trace.data, pen={'color': color, 'width': 2})
            rec_plot[t, c].setLabels(left=('Vm', 'V'))
            rec_plot[t, c].setLabels(bottom=('t', 's'))
            if deconv is True:
                rec_deconv = deconv_train(rec_pass_qc[:2])
                rec_deconv_grand = TSeriesList(rec_deconv[0]).mean()
                rec_ind_deconv_grand = TSeriesList(rec_deconv[1]).mean()
                #log_rec_plt.plot(rec_deconv_grand.time_values, rec_deconv_grand.data,
                #                    pen={'color': color, 'width': 2})
                rec_deconv_ind_grand2 = rec_ind_deconv_grand.copy(t0=delta + 0.2)
                #log_rec_plt.plot(rec_deconv_ind_grand2.time_values, rec_deconv_ind_grand2.data,
                #                    pen={'color': color, 'width': 2})
                [deconv_rec_plot[t, c].plot(ind.time_values, ind.data, pen=trace_color) for ind in rec_deconv[0]]
                [deconv_rec_plot[t, c].plot(rec.time_values, rec.data, pen=trace_color) for rec in rec_deconv[1]]
                deconv_rec_plot[t, c].plot(rec_ind_deconv_grand.time_values, rec_ind_deconv_grand.data,
                                    pen={'color': color, 'width': 2}, name=("n = %d" % n_synapses))
                deconv_rec_plot[t, c].plot(rec_deconv_grand.time_values, rec_deconv_grand.data,
                                    pen={'color': color, 'width': 2})
            summary_plot[c, 1].setLabels(left=('Norm Amp', ''))
            summary_plot[c, 1].setLabels(bottom=('Pulse Number', ''))
            f_color = pg.hsvColor(hue=hue, sat=float(t+0.5) / len(rec_t), val=1)
            summary_plot[c, 1].plot(rec_amp_grand/rec_amp_grand[0], name=('  %d ms' % delta), pen=f_color, symbol=symbols[t],
                                symbolBrush=f_color, symbolPen=None)
            if connection_types[c] not in rec_amp_summary.keys():
                rec_amp_summary[connection_types[c]] = []
                rec_uid[connection_types[c]] = []
            rec_amp_summary[connection_types[c]].append([delta, rec_amp])
            rec_uid[connection_types[c]].append([delta, rec_pass_qc[2]])

    # if sum(cache_change) > 0:
    #     write_cache(response_cache, cache_file)

print ('Exporting train pulse amplitudes and experiment IDs for further analysis')
write_cache([ind_amp_summary, rec_amp_summary], "train_amps_human.pkl")
write_cache([ind_uid, rec_uid], "expt_ids_human.pkl")