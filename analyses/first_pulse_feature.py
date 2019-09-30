import pyqtgraph as pg
import numpy as np
import csv
import sys
import argparse
from aisynphys.experiment_list import cached_experiments
from manuscript_figures import get_response, get_amplitude, response_filter, feature_anova, write_cache, trace_plot, \
    colors_human, colors_mouse, fail_rate, pulse_qc, feature_kw
from synapse_comparison import load_cache, summary_plot_pulse
from neuroanalysis.data import TSeriesList
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.fitting import fit_psp
from rep_connections import ee_connections, human_connections, no_include, all_connections, ie_connections, ii_connections, ei_connections
from aisynphys.synaptic_dynamics import DynamicsAnalyzer
from scipy import stats
import time
import pandas as pd

app = pg.mkQApp()
pg.dbg()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

# parser = argparse.ArgumentParser(description='Enter organism and type of connection you"d like to analyze ex: mouse ee (all mouse excitatory-'
#                 'excitatory). Alternatively enter a cre-type connection ex: sim1-sim1')
# parser.add_argument('--organism', dest='organism', help='Select mouse or human')
# parser.add_argument('--connection', dest='connection', help='Specify connections to analyze')
# args = vars(parser.parse_args(sys.argv[1:]))

all_expts = cached_experiments()
# manifest = {'Type': [], 'Connection': [], 'amp': [], 'latency': [],'rise':[], 'rise2080': [], 'rise1090': [], 'rise1080': [],
#             'decay': [], 'nrmse': [], 'CV': []}
# fit_qc = {'nrmse': 8, 'decay': 499e-3}

# if args['organism'] == 'mouse':
#     color_palette = colors_mouse
#     calcium = 'high'
#     age = '40-60'
#     sweep_threshold = 3
#     threshold = 0.03e-3
#     connection = args['connection']
#     if connection == 'ee':
#         connection_types = ee_connections.keys()
#     elif connection == 'ii':
#         connection_types = ii_connections.keys()
#     elif connection == 'ei':
#         connection_types = ei_connections.keys()
#     elif connection == 'ie':
#         connection_types = ie_connections.keys()
#     elif connection == 'all':
#         connection_types = all_connections.keys()
#     elif len(connection.split('-')) == 2:
#         c_type = connection.split('-')
#         if c_type[0] == '2/3':
#             pre_type = ('2/3', 'unknown')
#         else:
#             pre_type = (None, c_type[0])
#         if c_type[1] == '2/3':
#             post_type = ('2/3', 'unknown')
#         else:
#             post_type = (None, c_type[0])
#         connection_types = [(pre_type, post_type)]
# elif args['organism'] == 'human':
#     color_palette = colors_human
#     calcium = None
#     age = None
#     sweep_threshold = 5
#     threshold = None
#     connection = args['connection']
#     if connection == 'ee':
#         connection_types = human_connections.keys()
#     else:
#         c_type = connection.split('-')
#         connection_types = [((c_type[0], 'unknown'), (c_type[1], 'unknown'))]

plt = pg.plot()

# scale_offset = (-20, -20)
# scale_anchor = (0.4, 1)
# holding = [-65, -75]
# qc_plot = pg.plot()
# grand_response = {}
# expt_ids = {}
# feature_plot = None
# feature2_plot = PlotGrid()
# feature2_plot.set_shape(5,1)
# feature2_plot.show()
# feature3_plot = PlotGrid()
# feature3_plot.set_shape(1, 3)
# feature3_plot.show()
# amp_plot = pg.plot()
# synapse_plot = PlotGrid()
# synapse_plot.set_shape(len(connection_types), 1)
# synapse_plot.show()
# for c in range(len(connection_types)):
#     cre_type = (connection_types[c][0][1], connection_types[c][1][1])
#     target_layer = (connection_types[c][0][0], connection_types[c][1][0])
#     conn_type = connection_types[c]
#     #expt_list = all_expts.select(cre_type=cre_type, target_layer=target_layer, calcium=calcium, age=age)
#     color = color_palette[c]
#     grand_response[conn_type[0]] = {'trace': [], 'amp': [], 'latency': [], 'rise': [], 'dist': [], 'decay':[], 'CV': [], 'amp_measured': []}
#     expt_ids[conn_type[0]] = []
#     synapse_plot[c, 0].addLegend()
for expt in all_expts:
    if expt.connections is None:
        continue
    for pre, post in expt.connections:
        # cre_check = expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]
        # layer_check = expt.cells[pre].target_layer == target_layer[0] and expt.cells[post].target_layer == target_layer[1]
        # if cre_check is True and layer_check is True:
        pulse_response, artifact = get_response(expt, pre, post, analysis_type='pulse')
        stim_params = pulse_response.keys()
        pulse_traces = [sweep['response'] for sweep in pulse_response.values()]
        #     if threshold is not None and artifact > threshold:
        #         continue
        #     response_subset = response_filter(pulse_response, freq_range=[0, 50], holding_range=holding, pulse=True)
        #     if len(response_subset) >= sweep_threshold:
        #         qc_plot.clear()
        #         qc_list = pulse_qc(response_subset, baseline=1.5, pulse=None, plot=qc_plot)
        #         if len(qc_list) >= sweep_threshold:
        avg_trace, avg_amp, amp_sign, peak_t = get_amplitude(pulse_traces)
        #             if amp_sign is '-':
        #                 continue
        #             #print ('%s, %0.0f' %((expt.uid, pre, post), hold, ))
        #             all_amps = fail_rate(response_subset, '+', peak_t)
        #             cv = np.std(all_amps)/np.mean(all_amps)

                    # weight parts of the trace during fitting
        dt = avg_trace.dt
        weight = np.ones(len(avg_trace.data))*10.  #set everything to ten initially
        weight[int(10e-3/dt):int(12e-3/dt)] = 0.   #area around stim artifact
        weight[int(12e-3/dt):int(19e-3/dt)] = 30.  #area around steep PSP rise

        psp_fits = fit_psp(avg_trace,
                           xoffset=(14e-3, -float('inf'), float('inf')),
                           sign=amp_sign,
                           #                                           amp=avg_amp,
                           weight=weight)
        plt.clear()
        plt.plot(avg_trace.time_values, avg_trace.data, title=str([psp_fits.best_values['xoffset'], expt.uid, pre, post]))
        plt.plot(avg_trace.time_values, psp_fits.eval(), pen='g')
                    # avg_trace.t0 = -(psp_fits.best_values['xoffset'] - 10e-3)
                    # distance = expt.cells[pre].distance(expt.cells[post])
                    # grand_response[conn_type[0]]['CV'].append(cv)
                    # latency = psp_fits.best_values['xoffset'] - 10e-3
                    # rise = psp_fits.best_values['rise_time']
                    # decay = psp_fits.best_values['decay_tau']
                    # nrmse = psp_fits.nrmse()
                    # if nrmse < fit_qc['nrmse']:
                    #     grand_response[conn_type[0]]['latency'].append(psp_fits.best_values['xoffset'] - 10e-3)
                    #     max_x = np.argwhere(psp_fits.eval() == max(psp_fits.eval()))[0, 0]
                    #     min_x = int(psp_fits.best_values['xoffset'] / avg_trace.dt)
                    #     amp_20 = psp_fits.userkws['x'][np.argwhere(psp_fits.eval()[min_x:max_x] > psp_fits.best_values['amp']*0.2)[0,0]] + 10e-3
                    #     amp_10 = psp_fits.userkws['x'][np.argwhere(psp_fits.eval() > psp_fits.best_values['amp']*0.1)[0,0]]
                    #     amp_80 = psp_fits.userkws['x'][np.argwhere(psp_fits.eval()[min_x:max_x] < psp_fits.best_values['amp']*0.8)[-1,0]] +10e-3
                    #     amp_90 = psp_fits.userkws['x'][np.argwhere(psp_fits.eval()[: max_x] < psp_fits.best_values['amp']*0.9)[-1,0]]
                    #     rise2080 = amp_80-amp_20
                    #     rise1090 = amp_90-amp_10
                    #     rise1080 = amp_80-amp_10
                    #     grand_response[conn_type[0]]['rise'].append(rise2080)
                    #     manifest['rise2080'].append(rise2080)
                    #     manifest['rise1090'].append(rise1090)
                    #     manifest['rise1080'].append(rise1080)
                    #     manifest['rise'].append(rise)
                    # else:
                    #     grand_response[conn_type[0]]['latency'].append(np.nan)
                    #     grand_response[conn_type[0]]['rise'].append(np.nan)
                    #     manifest['rise2080'].append(None)
                    #     manifest['rise1090'].append(None)
                    #     manifest['rise1080'].append(None)
                    #     manifest['rise'].append(rise)
                    # grand_response[conn_type[0]]['trace'].append(avg_trace)
                    # grand_response[conn_type[0]]['amp'].append(psp_fits.best_values['amp'])
                    # grand_response[conn_type[0]]['amp_measured'].append(avg_amp)
                    # grand_response[conn_type[0]]['dist'].append(distance)
                    # expt_ids[conn_type[0]].append((pre, post, expt.uid, expt.source_id))
                    #
                    # manifest['Type'].append(cre_type)
                    # manifest['Connection'].append((expt.uid, pre, post))
                    # manifest['amp'].append(avg_amp)
                    # manifest['latency'].append(psp_fits.best_values['xoffset'] - 10e-3)
                    # manifest['nrmse'].append(psp_fits.nrmse())
                    # manifest['CV'].append(cv)
                    # manifest['decay'].append(psp_fits.best_values['decay_tau'])
                    #
                    # synapse_plot[c, 0].setTitle('First Pulse Response')
                    # if [expt.uid, pre, post] == all_connections[conn_type]:
                    #     trace_color = color + (30,)
                    # else:
                    #     trace_color = (0, 0, 0, 30)
                    # trace_plot(avg_trace, trace_color, plot=synapse_plot[c, 0], x_range=[0, 27e-3])
                    # app.processEvents()
#                    decay_response = response_filter(pulse_response, freq_range=[0, 20], holding_range=holding)
#                    qc_list = pulse_qc(response_subset, baseline=2, pulse=None, plot=qc_plot)
#                    if len(qc_list) >= sweep_threshold:
#                        avg_trace, avg_amp, amp_sign, peak_t = get_amplitude(qc_list)
#                        if amp_sign is '-':
#                            continue
#                        psp_fits = fit_psp(avg_trace, sign=amp_sign, yoffset=0, amp=avg_amp, method='leastsq', stacked = False,  fit_kws={})
#                        grand_response[conn_type[0]]['decay'].append(psp_fits.best_values['decay_tau'])
# if len(grand_response[conn_type[0]]['trace']) == 0:
#     continue
# if len(grand_response[conn_type[0]]['trace']) > 1:
#     grand_trace = TSeriesList(grand_response[conn_type[0]]['trace']).mean()
# else:
#     grand_trace = grand_response[conn_type[0]]['trace'][0]
# n_synapses = len(grand_response[conn_type[0]]['trace'])
# trace_plot(grand_trace, color={'color': color, 'width': 2}, plot=synapse_plot[c, 0], x_range=[0, 27e-3],
#            name=('%s, n = %d' % (connection_types[c], n_synapses)))
# synapse_plot[c, 0].hideAxis('bottom')
# # all_amps = np.hstack(np.asarray(grand_response[cre_type[0]]['fail_rate']))
# # y, x = np.histogram(all_amps, bins=np.linspace(0, 2e-3, 40))
# # synapse_plot[c, 1].plot(x, y, stepMode=True, fillLevel=0, brush='k')
# # synapse_plot[c, 1].setLabels(bottom=('Vm', 'V'))
# # synapse_plot[c, 1].setXRange(0, 2e-3)
# print ('%s kinetics n = %d' % (conn_type[0], len(grand_response[conn_type[0]]['latency'])))
# feature_list = (grand_response[conn_type[0]]['amp'], grand_response[conn_type[0]]['CV'], grand_response[conn_type[0]]['latency'],
#                 grand_response[conn_type[0]]['rise'])
# labels = (['Vm', 'V'], ['CV', ''], ['t', 's'], ['t', 's'])
# feature_plot = summary_plot_pulse(feature_list, labels, ('Amplitude', 'CV', 'Latency', 'Rise time'), c,
#                                   median=True, grand_trace=grand_trace, plot=feature_plot, color=color, name=connection_types[c])
# # feature2_plot[0, 0].plot(grand_response[conn_type[0]]['dist'], grand_response[conn_type[0]]['amp'], pen=None, symbol='o',
# #                          symbolSize=8, symbolBrush=color, symbolPeb='w')
# # feature2_plot[0, 0].setLabels(left=['mV', 'V'], bottom=['distance', 'm'])
# # feature2_plot[1, 0].plot(grand_response[conn_type[0]]['dist'], grand_response[conn_type[0]]['latency'], pen=None, symbol='o',
# #                          symbolSize=8, symbolBrush=color, symbolPeb='w')
# # feature2_plot[1, 0].setLabels(left=['ms', 's'], bottom=['distance', 'm'])
# # feature2_plot[2, 0].plot(grand_response[conn_type[0]]['dist'], grand_response[conn_type[0]]['rise'], pen=None, symbol='o',
# #                          symbolSize=8, symbolBrush=color, symbolPeb='w')
# # feature2_plot[2, 0].setLabels(left=['ms', 's'], bottom=['distance', 'm'])
# # feature2_plot[3, 0].plot(grand_response[conn_type[0]]['rise'], grand_response[conn_type[0]]['amp'], pen=None, symbol='o',
# #                          symbolSize=8, symbolBrush=color, symbolPeb='w')
# # feature2_plot[3, 0].setLabels(left=['amplitude', 'V'], bottom=['rise time', 's'])
# feature3_plot[0, 0].plot(grand_response[conn_type[0]]['latency'], grand_response[conn_type[0]]['amp'], pen=None, symbol='o',
#                          symbolSize=8, symbolBrush=color, symbolPen='w')
# feature3_plot[0, 0].setLabels(left=['Amp', 'V'], bottom=['Latency', 's'])
# feature3_plot[0, 1].plot(grand_response[conn_type[0]]['CV'], grand_response[conn_type[0]]['amp'], pen=None, symbol='o',
#                          symbolSize=8, symbolBrush=color, symbolPen='w')
# feature3_plot[0, 1].setLabels(left=['Amp', 'V'], bottom=['CV', ''])
# feature3_plot[0, 2].plot(grand_response[conn_type[0]]['CV'], grand_response[conn_type[0]]['latency'], pen=None, symbol='o',
#                          symbolSize=8, symbolBrush=color, symbolPen='w')
# feature3_plot[0, 2].setLabels(left=['Latency', 's'], bottom=['CV', ''])
# amp_plot.plot(grand_response[conn_type[0]]['amp'], grand_response[conn_type[0]]['amp_measured'], pen=None, symbol='o',
#              symbolSize=10, symbolBrush=color, symbolPen='w')
# amp_plot.setLabels(left=['Measured Amp', 'V'], bottom=['Fit Amp','V'])
# if c == len(connection_types) - 1:
#     x_scale = pg.ScaleBar(size=10e-3, suffix='s')
#     x_scale.setParentItem(synapse_plot[c, 0].vb)
#     x_scale.anchor(scale_anchor, scale_anchor, offset=scale_offset)
# amp_list = feature_anova('amp', grand_response)
# feature_kw('amp', grand_response)
# feature_kw('latency', grand_response)
# feature_kw('rise', grand_response)
# feature_kw('CV', grand_response)
# latency_list = feature_anova('latency', grand_response)
# rise_list = feature_anova('rise', grand_response)
# decay_list = feature_anova('decay', grand_response)

# if args['organism'] == 'human':
#     t, p = stats.ks_2samp(grand_response[('3', 'unknown')]['amp'], grand_response[('5', 'unknown')]['amp'])
#     print ('KS: Amp = %f' % p)
#     t, p = stats.ks_2samp(grand_response[('3', 'unknown')]['latency'],
#                            grand_response[('5', 'unknown')]['latency'])
#     print ('KS: Latency = %f' % p)
#     t, p = stats.ks_2samp(grand_response[('3', 'unknown')]['rise'],
#                            grand_response[('5', 'unknown')]['rise'])
#     print ('KS: Rise Time = %f' % p)
#     t, p = stats.ks_2samp(grand_response[('3', 'unknown')]['CV'],
#                           grand_response[('5', 'unknown')]['CV'])
#     print ('KS: CV = %f' % p)
# features = (amp_list, latency_list, rise_list, decay_list)
#
# #write_cache(expt_ids, 'pulse_expt_ids.pkl')
# #write_cache(features, 'pulse_features_human.pkl')
#
# df = pd.DataFrame(data=manifest)
# df = df[['Type', 'Connection', 'amp', 'latency', 'rise2080', 'rise1090', 'rise1080', 'decay', 'CV', 'nrmse']]
# writer = pd.ExcelWriter('Fig1_manifest_risefit.xlsx')
# df.to_excel(writer, 'Sheet1')
# writer.save()