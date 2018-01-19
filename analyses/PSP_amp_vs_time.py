"""
Question: how does PSP height change during the duration of the experiment.

"""


from __future__ import print_function, division

from collections import OrderedDict

import argparse
import sys
import pyqtgraph as pg
import os
import pickle
import pyqtgraph.multiprocess as mp
import numpy as np
import scipy.ndimage as ndi
import matplotlib.pyplot as plt

from multipatch_analysis.constants import INHIBITORY_CRE_TYPES
from multipatch_analysis.constants import EXCITATORY_CRE_TYPES
from multipatch_analysis.connection_detection import MultiPatchExperimentAnalyzer
from multipatch_analysis.synaptic_dynamics import DynamicsAnalyzer
from multipatch_analysis.experiment_list import cached_experiments
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import TraceList
from neuroanalysis.filter import bessel_filter
from neuroanalysis.event_detection import exp_deconvolve
from scipy import stats
import allensdk.core.json_utilities as ju
relative_path=os.path.dirname(os.getcwd())
sys.path.insert(1, os.path.join(relative_path))


def get_response(expt, pre, post, type='pulse'):
    '''This function was originally Stephanie's (analyses.manuscript_figures.get_response.
    I copy it here to play with it
    '''
    analyzer = DynamicsAnalyzer(expt, pre, post, method='deconv', align_to='spike')
    if type == 'pulse':
        response = analyzer.pulse_responses
        # pulse = 0  # only pull first pulse
        # response = {'data': [], 'dt': [], 'stim_param': []}
        # responses = analyzer.pulse_responses
        # for i,stim_params in enumerate(responses.keys()):
        #     resp = responses[stim_params]
        #     for trial in resp:
        #         r = trial[pulse]['response']
        #         r.meta['stim_params'] = stim_params
        #         response['data'].append(r.data)
        #         response['dt'].append(r.dt)
        #         response['stim_param'].append(r.meta['stim_params'])
    elif type == 'train':
        responses = analyzer.train_responses
        pulse_offset = analyzer.pulse_offsets
        response = {'responses': responses, 'pulse_offsets': pulse_offset}
    else:
        print ("Must select either pulse responses or train responses")
    if len(response) == 0:
        print ("No suitable data found for cell %d -> cell %d in expt %s" % (pre, post, expt.source_id))
        return response, None
    artifact = analyzer.cross_talk()
    return response, artifact

def measure_amp(trace, min_or_max, baseline=(6e-3, 8e-3), response=(13e-3, 17e-3)):
    '''get the max or min of the data in the trace object and subtract out the baseline
    at the specified times
    '''
    baseline = trace.time_slice(*baseline).data.mean()
    if min_or_max=='max':
        peak = trace.time_slice(*response).data.max()
    elif min_or_max=='min':
        peak = trace.time_slice(*response).data.min()
    else:
        raise Exception('Are you looking for min or max')
    return peak - baseline

def bin_data(time_list, data_list, bin_size=10):
    '''bins time series data in time bins with a specified size.
    Time must be bigger than 0.
    inputs
        time_list: list of arrays
            each array contains the times during the recording
        data_list: list of arrays
            data corresponding to time_list
        bin_size:
            specifies the size of the time bin
    returns:
        data_in_bins: list of numpy arrays
            each array corresponds to a time bin. Values in array are 
            values in the time bin
        time_bins: list
            values in list denote bin edges
        time_bins_middle: list
            values in list correspond to the center of time bins     
    '''
    max_time=max([max(tt) for tt in time_list])
    # make sure time and bin_size make sense
    if min([min(tt) for tt in time_list])<0: 
        raise Exception('time values should not be negative')
    if max_time<bin_size:
        raise Exception('bin size is bigger than max time')
    
    #specify the time bins 
    time_bins=np.arange(0, max_time+bin_size, bin_size) # this could potentially be broken depending on what the bin_size is
    
    if time_bins[-1] < max_time:
        raise Exception('Your largest time bin is less than max time.  Tweak your bin size.') 
    
    time_bin_middle=np.mean(np.array([np.append(time_bins, 0), np.append(0, time_bins)]), axis=0)[1:-1]  #time bin is defined by middle of bin
    data_in_bins=[np.array([]) for ii in range(len(time_bin_middle))] #initialize a data structure to receive data in bins
    
    #assign data to correct time bins
    for tt, ff in zip(time_list, data_list):
        assert len(tt)==len(ff)
        digits=np.digitize(tt, time_bins)-1 #note,digitize will assign values less than smallest timebin to a 0 index so -1 is used here
        for ii in range(len(digits)):
            data_in_bins[digits[ii]]=np.append(data_in_bins[digits[ii]], ff[ii])

    return data_in_bins, time_bins, time_bin_middle

def test_bin_data():
    '''meant to tests the bin_data() module but really just using a simple input so one can watch what is happening
    '''
    time_list=[np.array([1,2,3,4,5]), np.array([0, 3, 4.5]), np.array([1.1, 1.2, 2.3, 2.5, 2.8, 4, 4.3, 4.9])]
    data_list=[np.array([1,2,3,4,5]), np.array([6,7,8]), np.array([9,10,11,12,13,14,15, 16])]
    data, time_bins, time_mid_points=bin_data(time_list, data_list, bin_size=2.1)

def average_via_bins(time_list, data_list, bin_size=10):
    '''takes list of time arrays and corresponding list of data arrays and returns the average
    by placing the data in time bins
    '''
    data,time_bins, time_bin_middle=bin_data(time_list, data_list, bin_size)
    average_data=[]
    std_err_data=[]
    for bins in data:
        average_data.append(np.mean(bins))
        std_err_data.append(stats.sem(bins))
    assert len(average_data)==len(time_bin_middle), "data length doesn't match time length"
    return time_bin_middle, average_data, std_err_data
    

if __name__ == '__main__':
    app = pg.mkQApp()
    pg.dbg()
    
    # Load experiment index
    expts = cached_experiments()

    connection_list=[['rorb', 'rorb'],
                     ['tlx3', 'tlx3'],
                     ['ntsr1', 'ntsr1'],
                     ['L23pyr', 'L23pyr'],
                     ['sim1','sim1']]

    dictionary={}
    for synapic_pairs in connection_list:
        print(synapic_pairs)
    
        synapses = []
        for connection in expts.connection_summary():
            cells = connection['cells']
            expt = connection['expt']
            pre_synaptic=synapic_pairs[0]
            post_synaptic=synapic_pairs[1]
            if cells[0].cre_type == pre_synaptic and cells[1].cre_type == post_synaptic:
                synapses.append((expt, cells[0].cell_id, cells[1].cell_id))
        
        title_str= pre_synaptic+' to '+post_synaptic
        p = pg.plot(labels={'left': 'peak of synaptic deflection (mV)', 
                            'bottom': 'time since first recorded synapse (s)', 
                            'top':(title_str+' connections: progression of synaptic defection over an experiment')})    
        a = pg.plot(labels={'top':('average base-line subtracted first pulse synaptic deflection ('+ title_str+ ')'), 
                          'bottom': 'time (s)', 
                          'left':'voltage (mV)'}) 
        
        filtered=[]
        time_list=[]
        num_of_synapses=len(synapses)
        for i,syn in enumerate(synapses):
            expt, pre_id, post_id = syn
            analyzer = DynamicsAnalyzer(expt, pre_id, post_id, align_to='spike')
            
            # collect all first pulse responses
            responses = analyzer.amp_group
            if len(responses) == 0:
                print("Skipping %s %d %d; no responses" % (expt.uid, pre_id, post_id))
                continue
    
    #        plt.figure()
    #        for trace in responses.responses:
    #            plt.plot(trace.time_values, trace.data)
    #            plt.title('responses')
    #        plt.figure()
    #        for trace in responses.baselines:
    #            plt.plot(trace.time_values, trace.data)
    #            plt.title('baselines')         
    #        plt.show(block=False)
    
            # figure out whether the trough or peak of the average synaptic trace is bigger and if that corresponds to the synapse type.  
            # i.e. if it is an excitatory synapse we would expect the max defection to be positive
            average = responses.bsub_mean() #returns average synaptic response with the average baseline subtracted
            a.plot(average.time_values, average.data) #plot average of first pulse in each epoch of spikes of individual synapses
            max_peak = measure_amp(average, min_or_max='max', baseline=(6e-3, 8e-3), response=(12e-3, 16e-3))
            min_peak = measure_amp(average, min_or_max='min', baseline=(6e-3, 8e-3), response=(12e-3, 16e-3))
            max_min = "max" if abs(max_peak)> abs(min_peak) else "min"  #find whether the peak or trough of the first pulse average is larger 
            wrong_synapse_type_flag = False
            
            #if the largest deflection is in the negative direction 
            if max_min == "min" and expt.cells[pre_id].cre_type in EXCITATORY_CRE_TYPES: 
                print ("Whoa this synapse looks like a inhibitory when cre line would say it should be excitatory!!!" )  
                wrong_synapse_type_flag = True
            if max_min == "max" and expt.cells[pre_id].cre_type in INHIBITORY_CRE_TYPES: 
                print ("Whoa this synapse looks like a excitatory when cre line would say it should be inhibitory!!!" )    
                wrong_synapse_type_flag = True  
            
            # find the peak or trough of every potential event and plot their amplitude over time of the experiment
            peak=[]
            base=[]
            time=[]
            ordered=sorted(responses.responses, key=lambda rr:rr.start_time) #order the traces by time during the experiment
            for rr in ordered:
                peak.append(measure_amp(rr, min_or_max=max_min, baseline=(6e-3, 8e-3), response=(12e-3, 16e-3)))
                base.append(measure_amp(rr, min_or_max=max_min, baseline=(0e-3, 2e-3), response=(6e-3, 10e-3)))
                time.append(rr.start_time)        

    #        for trace in responses.responses:
    #            plt.plot(trace.time_values, trace.data)
    #            plt.title('responses')
    #        plt.figure()
    #        for trace in responses.baselines:
    #            plt.plot(trace.time_values, trace.data)
    #            plt.title('baselines')         
    #        plt.show()
            mean_base=np.mean(base)
            time=np.array(time) - time[0]
            time_list.append(time)
            peak_minus_base_average=np.array(peak)-mean_base
            filtered.append(ndi.gaussian_filter(peak_minus_base_average, 2))    
            
            if wrong_synapse_type_flag == False:
                p.plot(time, filtered[-1], pen=(i, len(synapses)*1.3))
            elif wrong_synapse_type_flag == True:
                p.plot(time, filtered[-1], pen=pg.mkPen(color=(i, len(synapses)*1.3),style=pg.QtCore.Qt.DashLine))
                            
            app.processEvents()
        
        #because times of events aren't all at the same time, time binning is needed to get average time course
        time_points, avg_data, std_err=average_via_bins(time_list, filtered, bin_size=60)
        p.plot(time_points, avg_data, pen=pg.mkPen(color='w', width=5)) #plots average of the data
        
        dictionary[title_str]={'time_points': time_points, 'avg_data':avg_data, 'std_err':std_err, 'num_of_synapses':num_of_synapses}
    ju.write("PSP_vs_time_output_data/psp_vs_time.json", dictionary)

    for key in dictionary.keys():
        plt.errorbar(dictionary[key]['time_points'], dictionary[key]['avg_data'],  yerr=dictionary[key]['std_err'], label=key)
    plt.title('average base-line subtracted first pulse synaptic deflection')
    plt.legend()
    plt.ylabel('voltage (mV)')
    plt.xlabel('time since first recorded synapse (s)')
    plt.show(block=False)

#        app.processEvents()    
#        pg.QtGui.QApplication.exec_()  
    
    plt.show()
        
    
  
    