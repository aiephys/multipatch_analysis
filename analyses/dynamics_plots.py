import pyqtgraph as pg
import numpy as np
import colorsys
from multipatch_analysis.experiment_list import cached_experiments
from manuscript_figures import cache_response, train_amp, induction_summary, recovery_summary, write_cache, get_response, train_qc
from synapse_comparison import load_cache
from neuroanalysis.data import TraceList
from neuroanalysis.ui.plot_grid import PlotGrid
from multipatch_analysis.synaptic_dynamics import RawDynamicsAnalyzer
from rep_connections import ee_connections, human_connections
pg.dbg()
app = pg.mkQApp()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

all_expts = cached_experiments()

connections = ee_connections[('sim1', 'sim1')]
connection_types = [('sim1', 'sim1')] #connections.keys()

calcium = 'high'
age = '40-60'
holding = [-68, -72]
freqs = [10, 20, 50, 100]
rec_t = [250, 500, 1000, 2000, 4000]
sweep_threshold = 5
qc_params = ('+', 0.5e-3)

# cache_file = 'train_response_cache.pkl'
# response_cache = load_cache(cache_file)
# cache_change = []
qc_plot = pg.plot()
ind_plot = PlotGrid()
ind_plot.set_shape(4, len(connection_types))
ind_plot.show()
rec_plot = PlotGrid()
rec_plot.set_shape(5,len(connection_types))
rec_plot.show()
summary_plot = PlotGrid()
summary_plot.set_shape(len(connection_types), 2)
summary_plot.show()
symbols = ['o', 's', 'd', '+', 't']
trace_color = (0, 0, 0, 5)

for c in range(len(connection_types)):
    induction_grand = {}
    recovery_grand = {}
    pulse_offset_ind = {}
    pulse_offset_rec = {}
    cre_type = connection_types[c]
    expt_list = all_expts.select(cre_type=cre_type, calcium=calcium, age=age)
    color = (c, len(connection_types)*1.3)
    summary_plot[c, 0].addLegend()
    summary_plot[c, 1].addLegend()
    for expt in expt_list:
        for pre, post in expt.connections:
            if expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]:
                train_response = get_response(expt, pre, post, type='train')
                induction_grand, pulse_offset_ind = induction_summary(train_response, freqs, holding, thresh=sweep_threshold,
                                                                ind_dict=induction_grand, offset_dict=pulse_offset_ind,
                                                                      qc=None)
                recovery_grand, pulse_offset_rec = recovery_summary(train_response, rec_t, holding, thresh=sweep_threshold,
                                                                rec_dict=recovery_grand, offset_dict=pulse_offset_rec)

    for f, freq in enumerate(freqs):
        if freq not in induction_grand.keys():
            print ("%d Hz not represented in data set" % freq)
            continue
        ind_offsets = pulse_offset_ind[freq]
        qc_plot.clear()
        ind_pass_qc = train_qc(induction_grand[freq], ind_offsets, amp=qc_params[1], sign=qc_params[0], plot=qc_plot)
        n_synapses = len(ind_pass_qc[0])
        if n_synapses > 0:
            induction_grand_trace = TraceList(ind_pass_qc[0]).mean()
            ind_rec_grand_trace = TraceList(ind_pass_qc[1]).mean()
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
            summary_plot[c, 0].setLabels(left=('Norm Amp', ''))
            summary_plot[c, 0].setLabels(bottom=('Pulse Number', ''))
            f_color = pg.hsvColor(1, sat=float(f+0.5)/len(freqs), val=1)
            summary_plot[c, 0].plot(ind_amp_grand/ind_amp_grand[0], name=('  %d Hz' % freq), pen=f_color, symbol=symbols[f],
                                symbolBrush=f_color, symbolPen=None)

    for t, delta in enumerate(rec_t):
        if delta not in recovery_grand.keys():
            print ("%d ms not represented in data set" % delta)
            continue
        recovery_grand_trace = TraceList(recovery_grand[delta][1]).mean()
        rec_ind_grand_trace = TraceList(recovery_grand[delta][0]).mean()
        n_synapses = len(recovery_grand[delta][1])
        rec_offsets = pulse_offset_rec[delta]
        rec_amp = train_amp(recovery_grand[delta], rec_offsets, '+')
        rec_amp_grand = np.mean(rec_amp, 0)
        if t == 0:
            rec_plot[t, c].setTitle(connection_types[c])
        if c == 0:
            label = pg.LabelItem('%d ms Recovery' % delta)
            label.setParentItem(rec_plot[t, c].vb)
            label.setPos(50, 0)
            summary_plot[c, 1].setTitle('Recovery')
        rec_plot[t, c].addLegend()
        [rec_plot[t, c].plot(ind.time_values, ind.data, pen=trace_color) for ind in recovery_grand[delta][0]]
        [rec_plot[t, c].plot(rec.time_values, rec.data, pen=trace_color) for rec in recovery_grand[delta][1]]
        rec_plot[t, c].plot(rec_ind_grand_trace.time_values, rec_ind_grand_trace.data, pen={'color': color, 'width': 2},
                            name=("n = %d" % n_synapses))
        rec_plot[t, c].plot(recovery_grand_trace.time_values, recovery_grand_trace.data, pen={'color': color, 'width': 2})
        rec_plot[t, c].setLabels(left=('Vm', 'V'))
        rec_plot[t, c].setLabels(bottom=('t', 's'))
        summary_plot[c, 1].setLabels(left=('Norm Amp', ''))
        summary_plot[c, 1].setLabels(bottom=('Pulse Number', ''))
        f_color = pg.hsvColor(color[0], sat=float(t+0.5) / len(rec_t), val=1)
        summary_plot[c, 1].plot(rec_amp_grand/rec_amp_grand[0], name=('  %d ms' % delta), pen=f_color, symbol=symbols[t],
                                symbolBrush=f_color, symbolPen=None)

    # if sum(cache_change) > 0:
    #     write_cache(response_cache, cache_file)