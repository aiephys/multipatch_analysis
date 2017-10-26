import pyqtgraph as pg
import numpy as np
from experiment_list import ExperimentList
from manuscript_figures import cache_response, get_amplitude, response_filter, feature_anova
from synapse_comparison import load_cache, summary_plot_pulse
from neuroanalysis.data import TraceList
from neuroanalysis.ui.plot_grid import PlotGrid
from connection_detection import fit_psp

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
age = '40-50'

cache_file = 'pulse_response_cache.pkl'
response_cache = load_cache(cache_file)

grand_response = {}
feature_plot = None
grid = PlotGrid()
grid.set_shape(4, 1)
synapse_plot = (grid[0, 0], grid[1, 0], grid[2, 0], grid[3, 0])
synapse_plot[0].grid = grid
grid.show()
for c in range(len(connection_types)):
    cre_type = connection_types[c].split('-')
    expt_list = all_expts.select(cre_type=cre_type, calcium=calcium, age=age)
    color = (c, len(connection_types)*1.3)
    grand_response[cre_type[0]] = {'trace': [], 'amp': [], 'latency': [], 'rise': [], 'decay':[]}
    synapse_plot[c].addLegend()
    for expt in expt_list:
        for pre, post in expt.connections:
            if expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]:
                pulse_response = cache_response(expt, pre, post, cache_file, response_cache, type='pulse')
                response_subset = response_filter(pulse_response, freq_range=[0, 50], holding_range=[-68, -72], pulse=True)
                if len(response_subset) >= 10:
                    avg_trace, avg_amp, amp_sign, peak_t = get_amplitude(response_subset)
                    if amp_sign is '-':
                        continue
                    psp_fits = fit_psp(avg_trace, sign=amp_sign, yoffset=0, amp=avg_amp, method='leastsq', fit_kws={})
                    avg_trace.t0 = -(psp_fits.best_values['xoffset'] - 10e-3)
                    grand_response[cre_type[0]]['latency'].append(psp_fits.best_values['xoffset'] - 10e-3)
                    grand_response[cre_type[0]]['rise'].append(psp_fits.best_values['rise_time'])
                    grand_response[cre_type[0]]['trace'].append(avg_trace)
                    grand_response[cre_type[0]]['amp'].append(avg_amp)
                    synapse_plot[c].setTitle('First Pulse Response')
                    synapse_plot[c].setLabels(left=('Vm', 'V'))
                    synapse_plot[c].setLabels(bottom=('t', 's'))
                    if [expt.uid, pre, post] == connections[cre_type[0]]:
                        trace_color = (255, 0, 255, 30)
                    else:
                        trace_color = (0, 0, 0, 30)
                    synapse_plot[c].plot(avg_trace.time_values, avg_trace.data, pen=trace_color)
                    synapse_plot[c].setXRange(0, 27e-3)
                    app.processEvents()
                decay_response = response_filter(pulse_response, freq_range=[0, 20], holding_range=[-68, -72])
                if len(response_subset) >= 10:
                    avg_trace, avg_amp, amp_sign, peak_t = get_amplitude(response_subset)
                    if amp_sign is '-':
                        continue
                    psp_fits = fit_psp(avg_trace, sign=amp_sign, yoffset=0, amp=avg_amp, method='leastsq', fit_kws={})
                    avg_trace.t0 = -(psp_fits.best_values['xoffset'] - 10e-3)
                    grand_response[cre_type[0]]['decay'].append(psp_fits.best_values['decay_tau'])
    grand_trace = TraceList(grand_response[cre_type[0]]['trace']).mean()
    grand_trace.t0 = 0
    n_synapses = len(grand_response[cre_type[0]]['trace'])
    synapse_plot[c].plot(grand_trace.time_values, grand_trace.data, pen={'color': color, 'width': 2},
                      name=('%s, n = %d' % (connection_types[c], n_synapses)))
    feature_list = (grand_response[cre_type[0]]['amp'], grand_response[cre_type[0]]['latency'], grand_response[cre_type[0]]['rise'])
    grand_amp = np.mean(np.array(grand_response[cre_type[0]]['amp']))
    grand_latency = np.mean(np.array(grand_response[cre_type[0]]['latency']))
    grand_rise = np.mean(np.array(grand_response[cre_type[0]]['rise']))
    grand_decay = np.mean(np.array(grand_response[cre_type[0]]['decay']))
    labels = (['Vm', 'V'], ['t', 's'], ['t', 's'])
    feature_plot = summary_plot_pulse(grand_trace, feature_list,(grand_amp, grand_latency, grand_rise), labels,
                                  ('Amplitude', 'Latency', 'Rise time'), c, plot=feature_plot,
                                  color=color, name=connection_types[c])

feature_anova('amp', grand_response)
feature_anova('latency', grand_response)
feature_anova('rise', grand_response)