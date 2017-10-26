import pyqtgraph as pg
import numpy as np
import colorsys
from experiment_list import ExperimentList
from manuscript_figures import cache_response, trace_avg, response_filter, train_amp
from synapse_comparison import load_cache
from neuroanalysis.data import TraceList
from neuroanalysis.ui.plot_grid import PlotGrid
from synaptic_dynamics import RawDynamicsAnalyzer
pg.dbg()
app = pg.mkQApp()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

all_expts = ExperimentList(cache='expts_cache.pkl')

connection_types = ['L23pyr-L23pyr']#, 'rorb-rorb', 'sim1-sim1', 'tlx3-tlx3']

connections = {'L23pyr': ['1501101571.17', 1, 5],
               'rorb': ['1502301827.80', 6, 8],
               'sim1': ['1490642434.41', 3, 5],
               'tlx3': ['1492545925.15', 8, 5]}

calcium = 'high'
age = '40-50'
holding = [-68, -72]
freqs = [10, 20, 50, 100]
rec_t = [250, 500, 1000, 2000, 4000]
sweep_threshold = 5

cache_file = 'train_response_cache.pkl'
response_cache = load_cache(cache_file)

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

for c in range(len(connection_types)):
    induction_grand = {}
    recovery_grand = {}
    pulse_offset_ind = {}
    pulse_offset_rec = {}
    cre_type = connection_types[c].split('-')
    expt_list = all_expts.select(cre_type=cre_type, calcium=calcium, age=age)
    color = (c, len(connection_types)*1.3)
    for expt in expt_list:
        for pre, post in expt.connections:
            if expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]:
                train_response = cache_response(expt, pre, post, cache_file, response_cache, type='train')
                if [expt.uid, pre, post] == connections[cre_type[0]]:
                    trace_color = (255, 0, 255, 10)
                else:
                    trace_color = (0, 0, 0, 5)
                for f, freq in enumerate(freqs):
                    induction_traces = {}
                    induction_traces['responses'] = response_filter(train_response['responses'], freq_range=[freq, freq],
                                                                    holding_range=holding, train=0)
                    induction_traces['pulse_offsets'] = response_filter(train_response['pulse_offsets'], freq_range=[freq, freq])
                    ind_rec_traces = response_filter(train_response['responses'], freq_range=[freq, freq], holding_range=holding,
                                                     train=1, delta_t=rec_t[0])
                    if len(induction_traces['responses']) >= sweep_threshold:
                        induction_avg = trace_avg(induction_traces['responses'])
                        ind_rec_avg = trace_avg(ind_rec_traces)
                        ind_rec_avg.t0 = induction_avg.time_values[-1] + 0.1
                        if freq not in induction_grand.keys():
                            induction_grand[freq] = [[], []]
                        induction_grand[freq][0].append(induction_avg)
                        induction_grand[freq][1].append(ind_rec_avg)
                        pulse_offset_ind[freq] = induction_traces['pulse_offsets']
                        ind_plot[f, c].plot(induction_avg.time_values, induction_avg.data, pen=trace_color)
                        ind_plot[f, c].plot(ind_rec_avg.time_values, ind_rec_avg.data, pen=trace_color)
                        app.processEvents()
                rec_ind_traces = response_filter(train_response['responses'], freq_range=[50, 50], holding_range=holding, train=0)
                for t, delta in enumerate(rec_t):
                    recovery_traces = {}
                    recovery_traces['responses'] = response_filter(train_response['responses'], freq_range=[50, 50],
                                                                   holding_range=holding, train=1, delta_t=delta)
                    recovery_traces['pulse_offsets'] = response_filter(train_response['pulse_offsets'], freq_range=[50, 50],
                                                      delta_t=delta)
                    if len(recovery_traces['responses']) >= sweep_threshold:
                        recovery_avg = trace_avg(recovery_traces['responses'])
                        rec_ind_avg = trace_avg(rec_ind_traces)
                        recovery_avg.t0 = (rec_ind_avg.time_values[-1]) + 0.1
                        if delta not in recovery_grand.keys():
                            recovery_grand[delta] = [[], []]
                        recovery_grand[delta][0].append(rec_ind_avg)
                        recovery_grand[delta][1].append(recovery_avg)
                        pulse_offset_rec[delta] = recovery_traces['pulse_offsets']
                        rec_plot[t, c].plot(rec_ind_avg.time_values, rec_ind_avg.data, pen=trace_color)
                        rec_plot[t, c].plot(recovery_avg.time_values, recovery_avg.data, pen=trace_color)
                        app.processEvents()

    for f, freq in enumerate(freqs):
        induction_grand_trace = TraceList(induction_grand[freq][0]).mean()
        ind_rec_grand_trace = TraceList(induction_grand[freq][1]).mean()
        n_synapses = len(induction_grand[freq][0])
        ind_offsets = pulse_offset_ind[freq]
        ind_amp = train_amp(induction_grand[freq], ind_offsets, '+')
        ind_amp_grand = np.mean(ind_amp, 0)
        if f == 0:
            ind_plot[f, c].setTitle(connection_types[c].split('-')[0])
        if c == 0:
            label = pg.LabelItem('%d Hz Induction' % freq)
            label.setParentItem(ind_plot[f, c].vb)
            label.setPos(50, 0)
            summary_plot[c, 0].setTitle('Induction')
        ind_plot[f, c].addLegend()
        ind_plot[f, c].plot(induction_grand_trace.time_values, induction_grand_trace.data, pen={'color': color, 'width': 2},
                            name=("n = %d" % n_synapses))
        ind_plot[f, c].plot(ind_rec_grand_trace.time_values, ind_rec_grand_trace.data, pen={'color': color, 'width': 2})
        ind_plot[f, c].setLabels(left=('Vm', 'V'))
        ind_plot[f, c].setLabels(bottom=('t', 's'))
        summary_plot[c, 0].setLabels(left=('Norm Amp', ''))
        summary_plot[c, 0].setLabels(bottom=('Pulse Number', ''))
        # f_color = colorsys.rgb_to_hsv(color[0], color[1], color[2])
        # summary_plot[c, 0].addLegend()
        # summary_plot[c, 0].plot(ind_amp_grand/ind_amp_grand[0], name=connection_types[c], symbol=symbols[f],
        #                        symbolPen=None, symbolBrush=(f_color, f/len(freqs)), pen=(f_color, f/len(freqs)))

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
        rec_plot[t, c].plot(rec_ind_grand_trace.time_values, rec_ind_grand_trace.data, pen={'color': color, 'width': 2},
                            name=("n = %d" % n_synapses))
        rec_plot[t, c].plot(recovery_grand_trace.time_values, recovery_grand_trace.data, pen={'color': color, 'width': 2})
        rec_plot[t, c].setLabels(left=('Vm', 'V'))
        rec_plot[t, c].setLabels(bottom=('t', 's'))
        summary_plot[c, 1].setLabels(left=('Norm Amp', ''))
        summary_plot[c, 1].setLabels(bottom=('Pulse Number', ''))
        summary_plot[c, 1].plot(rec_amp_grand/rec_amp_grand[0], symbol=symbols[t])
