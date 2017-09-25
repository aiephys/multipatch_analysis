import pyqtgraph as pg
import numpy as np
from experiment_list import ExperimentList
from manuscript_figures import cache_response, get_amplitude
from synapse_comparison import load_cache, summary_plot_pulse
from neuroanalysis.data import TraceList
from neuroanalysis.ui.plot_grid import PlotGrid

app = pg.mkQApp()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

all_expts = ExperimentList(cache='expts_cache.pkl')

connection_types = ['L23pyr-L23pyr', 'rorb-rorb', 'sim1-sim1', 'tlx3-tlx3']
calcium = 'high'
age = '40-50'

cache_file = 'pulse_response_cache.pkl'
response_cache = load_cache(cache_file)

grand_amp = {}
grand_kinetics = {}
amp_plot = None

for c in range(len(connection_types)):
    cre_type = connection_types[c].split('-')
    expt_list = all_expts.select(cre_type=cre_type, calcium=calcium, age=age)
    grand_amp[cre_type[0]] = {'trace': [], 'amp': []}
    grand_kinetics[cre_type[0]] = {'trace': []}
    grid = PlotGrid()
    grid.set_shape(1, 2)
    synapse_plot = (grid[0, 0], grid[0, 1])
    synapse_plot[0].grid = grid
    synapse_plot[0].addLegend()
    synapse_plot[1].addLegend()
    grid.show()
    for expt in expt_list:
        for pre, post in expt.connections:
            if expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]:
                pulse_amp, pulse_kinetics = cache_response(expt, pre, post, cache_file, response_cache)
                if len(pulse_amp) >= 10:
                    avg_amp_trace, avg_amp, amp_sign = get_amplitude(pulse_amp)
                    if amp_sign is '-':
                        continue
                    grand_amp[cre_type[0]]['trace'].append(avg_amp_trace)
                    grand_amp[cre_type[0]]['amp'].append(avg_amp)
                    synapse_plot[0].setTitle('First Pulse Amplitude')
                    synapse_plot[0].setLabels(left=('Vm', 'V'))
                    synapse_plot[0].setLabels(bottom=('t', 's'))
                    synapse_plot[0].plot(avg_amp_trace.time_values, avg_amp_trace.data, pen=(0, 0, 0, 10))
                if len(pulse_kinetics) >= 5:
                    avg_kin_trace, _, _ = get_amplitude(pulse_kinetics)
                    grand_kinetics[cre_type[0]]['trace'].append(avg_kin_trace)
                    synapse_plot[1].setTitle('First Pulse Kinetics')
                    synapse_plot[1].setLabels(left=('Vm', 'V'))
                    synapse_plot[1].setLabels(bottom=('t', 's'))
                    synapse_plot[1].plot(avg_kin_trace.time_values, avg_kin_trace.data, pen=(0, 0, 0, 10))
    grand_amp_trace = TraceList(grand_amp[cre_type[0]]['trace']).mean()
    n_synapses_amp = len(grand_amp[cre_type[0]]['trace'])
    grand_kin_trace = TraceList(grand_kinetics[cre_type[0]]['trace']).mean()
    n_synapses_kin = len(grand_kinetics[cre_type[0]]['trace'])
    synapse_plot[0].plot(grand_amp_trace.time_values, grand_amp_trace.data, name='%s, n = %d' % (connection_types[c], n_synapses_amp),
                      pen={'color': 'k', 'width': 2})
    synapse_plot[1].plot(grand_kin_trace.time_values, grand_kin_trace.data,
                      name='%s, n = %d' % (connection_types[c], n_synapses_kin), pen={'color': 'k', 'width': 2})
    grand_amp_scatter = np.mean(np.array(grand_amp[cre_type[0]]['amp']))
    amp_plot = summary_plot_pulse(grand_amp_trace, grand_amp[cre_type[0]]['amp'], grand_amp_scatter, c, plot=amp_plot,
                                  color=(c, len(connection_types)*1.3), name=connection_types[c])