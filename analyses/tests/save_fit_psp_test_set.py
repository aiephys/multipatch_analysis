"""Create a test set of data for testing the fit_psp function.  Uses Steph's 
original first_puls_feature.py code to filter out error causing data.

Example run statement
python save save_fit_psp_test_set.py --organism mouse --connection ee

Comment in the code that does the saving at the bottom

NOTE: the way this code runs fit_psp is currently out of date.  It will
need to be updated if use is desired.
"""


import pyqtgraph as pg
import numpy as np
import csv
import sys
import argparse
from multipatch_analysis.experiment_list import cached_experiments
from manuscript_figures import get_response, get_amplitude, response_filter, feature_anova, write_cache, trace_plot, \
    colors_human, colors_mouse, fail_rate, pulse_qc, feature_kw
from synapse_comparison import load_cache, summary_plot_pulse
from neuroanalysis.data import TraceList, Trace
from neuroanalysis.ui.plot_grid import PlotGrid
from multipatch_analysis.connection_detection import fit_psp
from rep_connections import ee_connections, human_connections, no_include, all_connections, ie_connections, ii_connections, ei_connections
from multipatch_analysis.synaptic_dynamics import DynamicsAnalyzer
from scipy import stats
import time
import pandas as pd
import json
import os

app = pg.mkQApp()
pg.dbg()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

parser = argparse.ArgumentParser(description='Enter organism and type of connection you"d like to analyze ex: mouse ee (all mouse excitatory-'
                'excitatory). Alternatively enter a cre-type connection ex: sim1-sim1')
parser.add_argument('--organism', dest='organism', help='Select mouse or human')
parser.add_argument('--connection', dest='connection', help='Specify connections to analyze')
args = vars(parser.parse_args(sys.argv[1:]))

all_expts = cached_experiments()
manifest = {'Type': [], 'Connection': [], 'amp': [], 'latency': [],'rise':[], 'rise2080': [], 'rise1090': [], 'rise1080': [],
            'decay': [], 'nrmse': [], 'CV': []}
fit_qc = {'nrmse': 8, 'decay': 499e-3}

if args['organism'] == 'mouse':
    color_palette = colors_mouse
    calcium = 'high'
    age = '40-60'
    sweep_threshold = 3
    threshold = 0.03e-3
    connection = args['connection']
    if connection == 'ee':
        connection_types = ee_connections.keys()
    elif connection == 'ii':
        connection_types = ii_connections.keys()
    elif connection == 'ei':
        connection_types = ei_connections.keys()
    elif connection == 'ie':
        connection_types == ie_connections.keys()
    elif connection == 'all':
        connection_types = all_connections.keys()
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
elif args['organism'] == 'human':
    color_palette = colors_human
    calcium = None
    age = None
    sweep_threshold = 5
    threshold = None
    connection = args['connection']
    if connection == 'ee':
        connection_types = human_connections.keys()
    else:
        c_type = connection.split('-')
        connection_types = [((c_type[0], 'unknown'), (c_type[1], 'unknown'))]

plt = pg.plot()

scale_offset = (-20, -20)
scale_anchor = (0.4, 1)
holding = [-65, -75]
qc_plot = pg.plot()
grand_response = {}
expt_ids = {}
feature_plot = None
feature2_plot = PlotGrid()
feature2_plot.set_shape(5,1)
feature2_plot.show()
feature3_plot = PlotGrid()
feature3_plot.set_shape(1, 3)
feature3_plot.show()
amp_plot = pg.plot()
synapse_plot = PlotGrid()
synapse_plot.set_shape(len(connection_types), 1)
synapse_plot.show()
for c in range(len(connection_types)):
    cre_type = (connection_types[c][0][1], connection_types[c][1][1])
    target_layer = (connection_types[c][0][0], connection_types[c][1][0])
    conn_type = connection_types[c]
    expt_list = all_expts.select(cre_type=cre_type, target_layer=target_layer, calcium=calcium, age=age)
    color = color_palette[c]
    grand_response[conn_type[0]] = {'trace': [], 'amp': [], 'latency': [], 'rise': [], 'dist': [], 'decay':[], 'CV': [], 'amp_measured': []}
    expt_ids[conn_type[0]] = []
    synapse_plot[c, 0].addLegend()
    for expt in expt_list:
        for pre, post in expt.connections:
            if [expt.uid, pre, post] in no_include:
                continue
            cre_check = expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]
            layer_check = expt.cells[pre].target_layer == target_layer[0] and expt.cells[post].target_layer == target_layer[1]
            if cre_check is True and layer_check is True:
                pulse_response, artifact = get_response(expt, pre, post, analysis_type='pulse')
                if threshold is not None and artifact > threshold:
                    continue
                response_subset, hold = response_filter(pulse_response, freq_range=[0, 50], holding_range=holding, pulse=True)
                if len(response_subset) >= sweep_threshold:
                    qc_plot.clear()
                    qc_list = pulse_qc(response_subset, baseline=1.5, pulse=None, plot=qc_plot)
                    if len(qc_list) >= sweep_threshold:
                        avg_trace, avg_amp, amp_sign, peak_t = get_amplitude(qc_list)
#                        if amp_sign is '-':
#                            continue
#                        #print ('%s, %0.0f' %((expt.uid, pre, post), hold, ))
#                        all_amps = fail_rate(response_subset, '+', peak_t)
#                        cv = np.std(all_amps)/np.mean(all_amps)
#                        
#                        # weight parts of the trace during fitting
                        dt = avg_trace.dt
                        weight = np.ones(len(avg_trace.data))*10.  #set everything to ten initially
                        weight[int(10e-3/dt):int(12e-3/dt)] = 0.   #area around stim artifact
                        weight[int(12e-3/dt):int(19e-3/dt)] = 30.  #area around steep PSP rise 
                        
                        # check if the test data dir is there and if not create it
                        test_data_dir='test_psp_fit'
                        if not os.path.isdir(test_data_dir):
                            os.mkdir(test_data_dir)
                            
                        save_dict={}
                        save_dict['input']={'data': avg_trace.data.tolist(),
                                            'dtype': str(avg_trace.data.dtype),
                                            'dt': float(avg_trace.dt),
                                            'amp_sign': amp_sign,
                                            'yoffset': 0, 
                                            'xoffset': 14e-3, 
                                            'avg_amp': float(avg_amp),
                                            'method': 'leastsq', 
                                            'stacked': False, 
                                            'rise_time_mult_factor': 10., 
                                            'weight': weight.tolist()} 
                        
                        # need to remake trace because different output is created
                        avg_trace_simple=Trace(data=np.array(save_dict['input']['data']), dt=save_dict['input']['dt']) # create Trace object
                        
                        psp_fits_original = fit_psp(avg_trace, 
                                           sign=save_dict['input']['amp_sign'], 
                                           yoffset=save_dict['input']['yoffset'], 
                                           xoffset=save_dict['input']['xoffset'], 
                                           amp=save_dict['input']['avg_amp'],
                                           method=save_dict['input']['method'], 
                                           stacked=save_dict['input']['stacked'], 
                                           rise_time_mult_factor=save_dict['input']['rise_time_mult_factor'], 
                                           fit_kws={'weights': save_dict['input']['weight']})  

                        psp_fits_simple = fit_psp(avg_trace_simple, 
                                           sign=save_dict['input']['amp_sign'], 
                                           yoffset=save_dict['input']['yoffset'], 
                                           xoffset=save_dict['input']['xoffset'], 
                                           amp=save_dict['input']['avg_amp'],
                                           method=save_dict['input']['method'], 
                                           stacked=save_dict['input']['stacked'], 
                                           rise_time_mult_factor=save_dict['input']['rise_time_mult_factor'], 
                                           fit_kws={'weights': save_dict['input']['weight']})  
                        print expt.uid, pre, post    
                        if psp_fits_original.nrmse()!=psp_fits_simple.nrmse():     
                            print '  the nrmse values dont match'
                            print '\toriginal', psp_fits_original.nrmse()
                            print '\tsimple', psp_fits_simple.nrmse()


                        
#                        save_dict['out']={}
#                        save_dict['out']['best_values']=psp_fits.best_values     
#                        save_dict['out']['best_fit']=psp_fits.best_fit.tolist()
#                        save_dict['out']['nrmse']=float(psp_fits.nrmse())
#                        with open(os.path.join(test_data_dir,expt.uid + '_' + str(pre) + '_' + str(post)+'NOTstacked.json'), 'w') as out_file:
#                            json.dump(save_dict, out_file)
