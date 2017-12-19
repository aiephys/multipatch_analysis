import pyqtgraph as pg
import numpy as np
from experiment_list import ExperimentList
from manuscript_figures import get_response, get_amplitude, response_filter, feature_anova, write_cache, trace_plot, color_palette
from synapse_comparison import load_cache, summary_plot_pulse
from neuroanalysis.data import TraceList
from neuroanalysis.ui.plot_grid import PlotGrid
from connection_detection import fit_psp
from rep_connections import ee_connections, human_connections

app = pg.mkQApp()
pg.dbg()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

all_expts = ExperimentList(cache='expts_cache.pkl')

connections = human_connections
connection_types = connections.keys()

calcium = None
age = None
sweep_threshold = 5
scale_offset = (-20, -20)
scale_anchor = (0.4, 1)

grand_response = {}
feature_plot = None
grid = PlotGrid()
grid.set_shape(len(connection_types), 1)
synapse_plot = []
grid.show()
for c in range(len(connection_types)):
    cre_type = connection_types[c]
    expt_list = all_expts.select(cre_type=cre_type, calcium=calcium, age=age)
    color = color_palette[c]
    grand_response[cre_type[0]] = {'trace': [], 'amp': [], 'latency': [], 'rise': [], 'decay':[]}
    synapse_plot.append(grid[c, 0])
    synapse_plot[c].addLegend()
    for expt in expt_list:
        for pre, post in expt.connections:
            if expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]:
                pulse_response = get_response(expt, pre, post, type='pulse')
                response_subset = response_filter(pulse_response, freq_range=[0, 50], holding_range=[-68, -72], pulse=True)
                if len(response_subset) >= sweep_threshold:
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
                    if [expt.uid, pre, post] == connections[cre_type]:
                        trace_color = (255, 0, 255, 30)
                    else:
                        trace_color = (0, 0, 0, 30)
                    trace_plot(avg_trace, trace_color, plot=synapse_plot[c], x_range=[0, 27e-3])
                    app.processEvents()
                decay_response = response_filter(pulse_response, freq_range=[0, 20], holding_range=[-68, -72])
                if len(response_subset) >= sweep_threshold:
                    avg_trace, avg_amp, amp_sign, peak_t = get_amplitude(response_subset)
                    if amp_sign is '-':
                        continue
                    psp_fits = fit_psp(avg_trace, sign=amp_sign, yoffset=0, amp=avg_amp, method='leastsq', fit_kws={})
                    grand_response[cre_type[0]]['decay'].append(psp_fits.best_values['decay_tau'])
    if len(grand_response[cre_type[0]]['trace']) == 0:
        continue
    if len(grand_response[cre_type[0]]['trace']) > 1:
        grand_trace = TraceList(grand_response[cre_type[0]]['trace']).mean()
        grand_trace.t0 = 0
    else:
        grand_trace = grand_response[cre_type[0]]['trace'][0]
    n_synapses = len(grand_response[cre_type[0]]['trace'])
    trace_plot(grand_trace, color={'color': color, 'width': 2}, plot=synapse_plot[c], x_range=[0, 27e-3],
               name=('%s, n = %d' % (connection_types[c], n_synapses)))
    synapse_plot[c].hideAxis('bottom')
    feature_list = (grand_response[cre_type[0]]['amp'], grand_response[cre_type[0]]['latency'], grand_response[cre_type[0]]['rise'])
    grand_amp = np.mean(np.array(grand_response[cre_type[0]]['amp']))
    grand_latency = np.mean(np.array(grand_response[cre_type[0]]['latency']))
    grand_rise = np.mean(np.array(grand_response[cre_type[0]]['rise']))
    grand_decay = np.mean(np.array(grand_response[cre_type[0]]['decay']))
    labels = (['Vm', 'V'], ['t', 's'], ['t', 's'])
    feature_plot = summary_plot_pulse(grand_trace, feature_list,(grand_amp, grand_latency, grand_rise), labels,
                                  ('Amplitude', 'Latency', 'Rise time'), c, plot=feature_plot,
                                  color=color, name=connection_types[c])
    if c == len(connection_types) - 1:
        x_scale = pg.ScaleBar(size=10e-3, suffix='s')
        x_scale.setParentItem(synapse_plot[c].vb)
        x_scale.anchor(scale_anchor, scale_anchor, offset=scale_offset)
feature_anova('amp', grand_response)
feature_anova('latency', grand_response)
feature_anova('rise', grand_response)
