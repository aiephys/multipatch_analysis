"""
Question: how does PSP height change with sweep number.

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

from constants import INHIBITORY_CRE_TYPES
from constants import EXCITATORY_CRE_TYPES
from connection_detection import MultiPatchExperimentAnalyzer
from synaptic_dynamics import DynamicsAnalyzer
from experiment_list import ExperimentList
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import TraceList
from neuroanalysis.filter import bessel_filter
from neuroanalysis.event_detection import exp_deconvolve


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
    '''
    max_time=max([max(tt) for tt in time_list])
    if min([min(tt) for tt in time_list])<0: 
        raise Exception('time values should not be negative')
    if max_time<bin_size:
        raise Exception('bin size is bigger than max time')
    time_bins=np.arange(0, max_time+bin_size, bin_size) # this could potentially be broken depending on what the bin_size is
    if time_bins[-1] < max_time:
        raise Exception('Your largest time bin is less than max time.  Tweak your bin size.') 
    time_bin_middle=np.mean(np.array([np.append(time_bins, 0), np.append(0, time_bins)]), axis=0)[1:-1] 
    data_in_bins=[np.array([]) for ii in range(len(time_bin_middle))] #initialize a data structure to receive data in bins
    for tt, ff in zip(time_list, data_list):
        assert len(tt)==len(ff)
        digits=np.digitize(tt, time_bins)-1 #digitize will assign values less than smalled timebin to a 0 index so -1 is used here
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
    for bins in data:
        average_data.append(np.mean(bins))
    assert len(average_data)==len(time_bin_middle), "data length doesn't match time length"
    return time_bin_middle, average_data
    

if __name__ == '__main__':
    app = pg.mkQApp()
    
    # Load experiment index
    cache_file = 'expts_cache.pkl'
    expts = ExperimentList(cache=cache_file)

    synapses = []
    for connection in expts.connection_summary():
        cells = connection['cells']
        expt = connection['expt']
        if cells[0].cre_type == 'sim1' and cells[1].cre_type == 'sim1':
            synapses.append((expt, cells[0].cell_id, cells[1].cell_id))
            
    p = pg.plot(labels={'left': 'peak of synaptic deflection (mV)', 'bottom': 'time since first recorded synapse (s)', 'top':'sim1 to sim1 connections'})    
    a=pg.plot(labels={'top':'average base-line subtraced first pluse synaptic deflection (sim1 to sim1)', 'bottom': 'time (s)', 'left':'voltage (mV)'}) 

    filtered=[]
    time_list=[]
    for i,syn in enumerate(synapses):
        expt, pre_id, post_id = syn
        analyzer = DynamicsAnalyzer(expt, pre_id, post_id, align_to='spike')
        
        # collect all first pulse responses
        responses = analyzer.amp_group
        if len(responses) == 0:
            print("Skipping %s %d %d; no responses" % (expt.uid, pre_id, post_id))
            continue

        # figure out whether the trough or peak of the average synaptic trace is bigger and if that corresponds to the synapse type.  
        # i.e. if it is an excitatory synapse we would expect the max defection to be negative
        average = responses.bsub_mean() #returns average synaptic response with the average baseline subtracted
        a.plot(average.time_values, average.data) #plot average of first pulse of individual synapses
        max_peak = measure_amp(average, min_or_max='max', baseline=(6e-3, 8e-3), response=(12e-3, 16e-3))
        min_peak = measure_amp(average, min_or_max='min', baseline=(6e-3, 8e-3), response=(12e-3, 16e-3))
        max_min = "max" if abs(max_peak)> abs(min_peak) else "min"  #find whether the peak or trough is larger
        wrong_synapse_type_flag = False
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
    time_points, avg_data=average_via_bins(time_list, filtered, bin_size=60)
    p.plot(time_points, avg_data, pen=pg.mkPen(color='w', width=5)) #plots average of the data
    

    app.processEvents()    
    pg.QtGui.QApplication.exec_()    
    