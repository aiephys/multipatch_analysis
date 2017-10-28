import pyqtgraph as pg
import numpy as np
import colorsys
from experiment_list import ExperimentList
from manuscript_figures import cache_response, train_amp, induction_summary, recovery_summary, write_cache
from synapse_comparison import load_cache
from neuroanalysis.data import TraceList
from neuroanalysis.ui.plot_grid import PlotGrid
from synaptic_dynamics import RawDynamicsAnalyzer
pg.dbg()
app = pg.mkQApp()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

all_expts = ExperimentList(cache='expts_cache.pkl')

connection_types = ['L23pyr-L23pyr', 'rorb-rorb', 'sim1-sim1', 'tlx3-tlx3']

connections = {'L23pyr': ['1501101571.17', 1, 5],
               'rorb': ['1502301827.80', 6, 8],
               'sim1': ['1490642434.41', 3, 5],
               'tlx3': ['1492545925.15', 8, 5]}

calcium = 'high'
age = '40-60'
holding = [-68, -72]
freqs = [10, 20, 50, 100]
rec_t = [250, 500, 1000, 2000, 4000]
sweep_threshold = 5

cache_file = 'train_response_cache.pkl'
response_cache = load_cache(cache_file)
cache_change = []

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
    cre_type = connection_types[c].split('-')
    expt_list = all_expts.select(cre_type=cre_type, calcium=calcium, age=age)
    color = (c, len(connection_types)*1.3)
    summary_plot[c, 0].addLegend()
    summary_plot[c, 1].addLegend()
    for expt in expt_list:
        for pre, post in expt.connections:
            if expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]:
                train_response, change = cache_response(expt, pre, post, response_cache, type='train')
                cache_change.append(change)
                induction_grand, pulse_offset_ind = induction_summary(train_response, freqs, holding, thresh=sweep_threshold,
                                                                ind_dict=induction_grand, offset_dict=pulse_offset_ind)
                recovery_grand, pulse_offset_rec = recovery_summary(train_response, rec_t, holding, thresh=sweep_threshold,
                                                                rec_dict=recovery_grand, offset_dict=pulse_offset_rec)

    for f, freq in enumerate(freqs):
        induction_grand_trace = TraceList(induction_grand[freq][0]).mean()
        ind_rec_grand_trace = TraceList(induction_grand[freq][1]).mean()
        n_synapses = len(induction_grand[freq][0])
        ind_offsets = pulse_offset_ind[freq]
        ind_amp = train_amp(induction_grand[freq], ind_offsets, '+')
        ind_amp_grand = np.mean(ind_amp, 0)
        if f == 0:
            ind_plot[f, c].setTitle(connection_types[c].split('-')[0])
            type = pg.LabelItem('%s' % connection_types[c])
            type.setParentItem(summary_plot[c, 0])
            type.setPos(50, 0)
        if c == 0:
            label = pg.LabelItem('%d Hz Induction' % freq)
            label.setParentItem(ind_plot[f, c].vb)
            label.setPos(50, 0)
            summary_plot[c, 0].setTitle('Induction')
        ind_plot[f, c].addLegend()
        [ind_plot[f, c].plot(ind.time_values, ind.data, pen=trace_color) for ind in induction_grand[freq][0]]
        [ind_plot[f, c].plot(rec.time_values, rec.data, pen=trace_color) for rec in induction_grand[freq][1]]
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
        recovery_grand_trace = TraceList(recovery_grand[delta][1]).mean()
        rec_ind_grand_trace = TraceList(recovery_grand[delta][0]).mean()
        n_synapses = len(recovery_grand[delta][1])
        rec_offsets = pulse_offset_rec[delta]
        rec_amp = train_amp(recovery_grand[delta], rec_offsets, '+')
        rec_amp_grand = np.mean(rec_amp, 0)
        if t == 0:
            rec_plot[t, c].setTitle(connection_types[c].split('-')[0])
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

    if sum(cache_change) > 0:
        write_cache(response_cache, cache_file)