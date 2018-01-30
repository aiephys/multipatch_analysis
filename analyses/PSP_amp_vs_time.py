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

import datetime
from multipatch_analysis.constants import INHIBITORY_CRE_TYPES
from multipatch_analysis.constants import EXCITATORY_CRE_TYPES
from multipatch_analysis.connection_detection import MultiPatchExperimentAnalyzer
from multipatch_analysis.synaptic_dynamics import DynamicsAnalyzer
from multipatch_analysis.experiment_list import cached_experiments
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import TraceList, PatchClampRecording
from neuroanalysis.filter import bessel_filter
from neuroanalysis.event_detection import exp_deconvolve
from scipy import stats
import allensdk.core.json_utilities as ju
relative_path=os.path.dirname(os.getcwd())
sys.path.insert(1, os.path.join(relative_path))


def response_filter(response, freq_range=None, holding_range=None, pulse=False, train=None, delta_t=None):
    '''this is from Stephanie's manuscript_figures.py
    '''
    #plot = pg.plot()
    new_responses = []
    for stim_params, trials in response.items():
        ind_freq, rec_t, holding = stim_params
        holding = holding * 1e3
        rec_t = int(np.round(rec_t * 1e3, -1))
        if freq_range is not None and (ind_freq < freq_range[0] or ind_freq > freq_range[1]):
            continue
        if holding_range is not None and (holding > holding_range[0] or holding < holding_range[1]):
            continue
        if delta_t is not None and rec_t != delta_t:
            continue
        if pulse is True:
            for trial in trials:
                new_responses.append(trial[0]['response'])
        elif train is not None:
            for trial in trials[train].responses:
                new_responses.append(trial)
        else:
            new_responses.append(trials)
        # plot.plot(response[trial].time_values, response[trial].data)
        # app.processEvents()
    return new_responses

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

def bin_data(the_list, data_list, bin_size=10):
    '''bins time series data in time bins with a specified size.
    Time must be bigger than 0.
    inputs
        the_list: list of arrays
            each array contains the times during the recording
        data_list: list of arrays
            data corresponding to the_list
        bin_size:
            specifies the size of the time bin
    returns:
        data_in_bins: list of numpy arrays
            each array corresponds to a time bin. Values in array are 
            values in the time bin
        the_bins: list
            values in list denote bin edges
        middle_of_bins: list
            values in list correspond to the center of time bins     
    '''
    max_value=max([max(tt) for tt in the_list])
    # make sure time and bin_size make sense
    if min([min(tt) for tt in the_list])<0: 
        raise Exception('time values should not be negative')
    if max_value<bin_size:
        raise Exception('bin size is bigger than max time')
    
    #specify the time bins 
    the_bins=np.arange(0, max_value+bin_size, bin_size) # this could potentially be broken depending on what the bin_size is
    
    if the_bins[-1] < max_value:
        raise Exception('Your largest time bin is less than max time.  Tweak your bin size.') 
    
    middle_of_bins=np.mean(np.array([np.append(the_bins, 0), np.append(0, the_bins)]), axis=0)[1:-1]  #time bin is defined by middle of bin
    data_in_bins=[np.array([]) for ii in range(len(middle_of_bins))] #initialize a data structure to receive data in bins
    
    #assign data to correct time bins
    for tt, ff in zip(the_list, data_list):
        assert len(tt)==len(ff)
        digits=np.digitize(tt, the_bins)-1 #note,digitize will assign values less than smallest timebin to a 0 index so -1 is used here
        for ii in range(len(digits)-1):  #note I just added this -1 for the sweep indexing so not sure if it makes sense 
            print (ii, len(digits))
            print (data_in_bins[digits[ii]], ff[ii])
            data_in_bins[digits[ii]]=np.append(data_in_bins[digits[ii]], ff[ii])

    return data_in_bins, the_bins, middle_of_bins

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
    
def check_synapse(expt, cells):
    '''checks if a synapse meets the requirements and if so, it appends it to the synapse 
    dictionary.  A synapses list must be initialized before calling this function
    inputs:
        expt: object
            object obtained from cached_experiments.connection_summary[*]['expt']
        cells: object
            object obtained from cached_experiments.connection_summary[*]['cells']    
    output:
        returns nothing but appends info to the synapse list
    '''
    try: #needed because pyqt is breaking on datetime sometimes
        if expt.expt_info['solution']=='2mM Ca & Mg':
            synapses.append((expt, cells[0].cell_id, cells[1].cell_id))
    except:
        pass
        
    
if __name__ == '__main__':
    app = pg.mkQApp()
    pg.dbg()
    
    # Load experiment index
    expts = cached_experiments()
#    expts.select(calcium='high')  #this is throwing datetime errors

    connection_list=[
                     ['rorb', 'rorb'],
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
            if pre_synaptic=='L23pyr':
                if cells[0].target_layer=='2/3' and cells[1].target_layer=='2/3':
                    check_synapse(expt, cells)
            else:
                if cells[0].cre_type == pre_synaptic and cells[1].cre_type == post_synaptic:
                    check_synapse(expt, cells)
        
        title_str= pre_synaptic+' to '+post_synaptic
        time_vs_psp_plot = pg.plot(labels={'left': 'peak of synaptic deflection (V)', 
                            'bottom': 'time since first recorded synapse (s)', 
                            'top':(title_str+' connections: progression of synaptic defection over an experiment')})    
        ave_psp_plot = pg.plot(labels={'top':('average base-line subtracted first pulse synaptic deflection ('+ title_str+ ')'), 
                          'bottom': 'time (s)', 
                          'left':'voltage (V)'}) 
        sweep_vs_psp_plot = pg.plot(labels={'left': 'peak of synaptic deflection (V)', 
                            'bottom': 'sweep number', 
                            'top':(title_str+' connections: progression of synaptic defection over an experiment')})  
        
        raw=[]
        filtered=[]
        time_list=[]
        sweep_number_list=[]
        num_of_synapses=0
        for i,syn in enumerate(synapses):
            expt, pre_id, post_id = syn
            analyzer = DynamicsAnalyzer(expt, pre_id, post_id, align_to='spike')
            
            # collect all first pulse responses
            amp_responses = analyzer.amp_group
            if len(amp_responses) == 0:
                print("Skipping %s %d %d; no responses" % (expt.uid, pre_id, post_id))
                continue
            
            
            # this doesnt work yet it has a pyqt error
#            train_responses=analyzer.train_responses
#            good_holding_responses=response_filter(train_responses,  holding_range=[-68, -72], pulse=True)
#            
            
            
# some other options to potentially be able to choose from            
#            responses = analyzer.train_responses
#            pulse_offset = analyzer.pulse_offsets
#            response = analyzer.pulse_responses
    
    #        plt.figure()
    #        for trace in amp_responses.responses:
    #            plt.plot(trace.time_values, trace.data)
    #            plt.title('responses')
    #        plt.figure()
    #        for trace in amp_responses.baselines:
    #            plt.plot(trace.time_values, trace.data)
    #            plt.title('baselines')         
    #        plt.show(block=False)
    
            # figure out whether the trough or peak of the average synaptic trace is bigger and if that corresponds to the excitation of the neurons.  
            # i.e. if it is an excitatory synapse we would expect the max defection to be positive
            average = amp_responses.bsub_mean() #returns average synaptic response with the average baseline subtracted
            max_peak = measure_amp(average, min_or_max='max', baseline=(6e-3, 8e-3), response=(12e-3, 16e-3))
            min_peak = measure_amp(average, min_or_max='min', baseline=(6e-3, 8e-3), response=(12e-3, 16e-3))
            max_min = "max" if abs(max_peak)> abs(min_peak) else "min"  #find whether the peak or trough of the first pulse average is larger 
            wrong_synapse_type_flag = False
            if max_min == "min" and expt.cells[pre_id].cre_type in EXCITATORY_CRE_TYPES: 
                print ("Whoa this synapse looks inhibitory when cre line would say it should be excitatory!!!" )  
                wrong_synapse_type_flag = True
            if max_min == "max" and expt.cells[pre_id].cre_type in INHIBITORY_CRE_TYPES: 
                print ("Whoa this synapse looks excitatory when cre line would say it should be inhibitory!!!" )    
                wrong_synapse_type_flag = True  
       
            # find the peak or trough of every potential event and plot their amplitude over time of the experiment
            peak=[]
            base=[]
            time=[]
            holding_potential=[]
            sweep_number=[]
            ordered=sorted(amp_responses.responses, key=lambda rr:rr.start_time) #order the traces by time during the experiment
            for rr in ordered:
                peak.append(measure_amp(rr, min_or_max=max_min, baseline=(6e-3, 8e-3), response=(12e-3, 16e-3)))
                base.append(measure_amp(rr, min_or_max=max_min, baseline=(0e-3, 2e-3), response=(6e-3, 10e-3)))
                time.append(rr.start_time)  
                holding_potential.append(rr.parent.holding_potential)
                sweep_number.append(rr.parent.parent._sweep_id)
                print ('for each pulse of a synapse: peak', peak[-1], 'base', base[-1], 'time', time[-1], 'holding potential', holding_potential[-1], 'sweep_number', sweep_number[-1])      

    #        for trace in amp_responses.responses:
    #            plt.plot(trace.time_values, trace.data)
    #            plt.title('responses')
    #        plt.figure()
    #        for trace in amp_responses.baselines:
    #            plt.plot(trace.time_values, trace.data)
    #            plt.title('baselines')         
    #        plt.show()
    
            # check if holding potential is within a desired range
            holding=np.mean(holding_potential) # average holding potential across plots
            print ('holding potential is', holding)
            if holding>-0.072 and holding<-0.068:
                holding_good_flag=True
            else:
                holding_good_flag=False
            print ('\tholding potential flag set to ', holding_good_flag)

            # if a neuron passes the criterea for inclusion add the synapse to the list
            mean_base=np.mean(base) # average base across pulses of a synapse
            time=np.array(time) - time[0] # remap time basis to be in reference to start of experiment
            peak_minus_base_average=np.array(peak)-mean_base # take each peak and put it in reference to the average base
            smoothed=ndi.gaussian_filter(peak_minus_base_average, 2) # 
            if wrong_synapse_type_flag == False and holding_good_flag ==True:
                print('recording synapse')
                time_list.append(time)  
                filtered.append(smoothed)
                raw.append(peak_minus_base_average)
                sweep_number_list.append(sweep_number)
                ave_psp_plot.plot(average.time_values, average.data, pen=pg.mkPen(color=(0, 128, 0))) #plot average of first pulse in each epoch of spikes of individual synapses
                time_vs_psp_plot.plot(time, smoothed, pen=pg.mkPen(color=(0, 128, 0))) # (i, len(synapses)*1.3))
                sweep_vs_psp_plot.plot(sweep_number, smoothed, pen=pg.mkPen(color=(0, 128, 0)))
                num_of_synapses=num_of_synapses+1
            else: 
                print ('wrong_synapse_type_flag', wrong_synapse_type_flag)
                print ('holding_good_flag', holding_good_flag)
                if wrong_synapse_type_flag==True and holding_good_flag==False:
                    ave_psp_plot.plot(average.time_values, average.data, pen=pg.mkPen(color=(255, 0, 0))) #plot average of first pulse in each epoch of spikes of individual synapses
                    time_vs_psp_plot.plot(time, smoothed, pen=pg.mkPen(color=(255, 0, 0)))#pen=pg.mkPen(color=(i, len(synapses)*1.3),style=pg.QtCore.Qt.DashDotLine))             
                    sweep_vs_psp_plot.plot(sweep_number, smoothed, pen=pg.mkPen(color=(0, 128, 0)))
                elif wrong_synapse_type_flag==True and holding_good_flag==True:
                    ave_psp_plot.plot(average.time_values, average.data, pen=pg.mkPen(color=(0,191,255))) #plot average of first pulse in each epoch of spikes of individual synapses
                    time_vs_psp_plot.plot(time, smoothed, pen=pg.mkPen(color=(0,191,255)))#pen=pg.mkPen(color=(i, len(synapses)*1.3),style=pg.QtCore.Qt.DashDotLine))
                    sweep_vs_psp_plot.plot(sweep_number, smoothed, pen=pg.mkPen(color=(0, 191, 255)))
                elif wrong_synapse_type_flag==False and holding_good_flag==False:
                    ave_psp_plot.plot(average.time_values, average.data, pen=pg.mkPen(color=(138,43,226))) #plot average of first pulse in each epoch of spikes of individual synapses
                    time_vs_psp_plot.plot(time, smoothed, pen=pg.mkPen(color=(138,43,226)))#pen=pg.mkPen(color=(i, len(synapses)*1.3),style=pg.QtCore.Qt.DashDotLine))
                    sweep_vs_psp_plot.plot(sweep_number, smoothed, pen=pg.mkPen(color=(138, 43, 226)))
                else:
                    raise Exception('This flag combo doesnt exist')
            app.processEvents()
        
        #because times of events aren't all at the same time, time binning is needed to get average time course
        time_points, time_avg_data, time_std_err=average_via_bins(time_list, raw, bin_size=60)
        time_vs_psp_plot.plot(time_points, time_avg_data, pen=pg.mkPen(color='w', width=5)) #plots average of the data
        
        sweeps, sweep_avg_data, sweep_std_err=average_via_bins(sweep_number_list, raw, bin_size=5)
        sweep_vs_psp_plot.plot(sweeps, sweep_avg_data, pen=pg.mkPen(color='w', width=5)) #plots average of the data
        
        dictionary[title_str]={'time_points': time_points, 
                               'time_avg_data':time_avg_data, 
                               'time_std_err':time_std_err, 
                               'num_of_synapses':num_of_synapses,
                               'sweeps':sweeps,
                               'sweep_avg_data':sweep_avg_data,
                               'sweep_std_err':sweep_std_err}
        
    ju.write("PSP_vs_time_output_data/goodpsp_vs_time or_sweep_1_29_18.json", dictionary)

    plt.figure()
    for key in dictionary.keys():
        plt.errorbar(dictionary[key]['time_points'], dictionary[key]['time_avg_data'],  yerr=dictionary[key]['time_std_err'], label=key+', n='+str(dictionary[key]['num_of_synapses']))
    plt.title('average base-line subtracted first pulse synaptic deflection')
    plt.legend(loc=4)
    plt.ylabel('voltage (mV)')
    plt.xlabel('time since first recorded synapse (s)')

    plt.figure()
    for key in dictionary.keys():
        plt.errorbar(dictionary[key]['sweeps'], dictionary[key]['sweep_avg_data'],  yerr=dictionary[key]['sweep_std_err'], label=key+', n='+str(dictionary[key]['num_of_synapses']))
    plt.title('average base-line subtracted first pulse synaptic deflection')
    plt.legend(loc=4)
    plt.ylabel('voltage (mV)')
    plt.xlabel('sweep number')    
    
    
    plt.show(block=False)

#        app.processEvents()    
#        pg.QtGui.QApplication.exec_()  
    
    plt.show()
        
    
  
    