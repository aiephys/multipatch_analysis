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
import scipy.stats as stats

colors={}
colors['correct_amp_good_HP']=(0, 128, 0) #green
colors['correct_amp_bad_HP']=(138,43,226) #purple
colors['wrong_amp_good_HP']=(0,191,255) #cyan
colors['wrong_amp_bad_HP']=(255, 0, 0) #red

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
        individual_first_peaks = trace.time_slice(*response).data.max()
    elif min_or_max=='min':
        individual_first_peaks = trace.time_slice(*response).data.min()
    else:
        raise Exception('Are you looking for min or max?')
    return individual_first_peaks - baseline

def bin_data(the_list, data_list, bin_size=10):
    '''bins individual_start_times series data in individual_start_times bins with a specified size.
    Time must be bigger than 0.
    inputs
        the_list: list of arrays
            each array contains the times during the recording
        data_list: list of arrays
            data corresponding to the_list
        bin_size:
            specifies the size of the individual_start_times bin
    returns:
        data_in_bins: list of numpy arrays
            each array corresponds to a individual_start_times bin. Values in array are 
            values in the individual_start_times bin
        the_bins: list
            values in list denote bin edges
        middle_of_bins: list
            values in list correspond to the center of individual_start_times bins     
    '''
    max_value=max([max(tt) for tt in the_list])
    # make sure individual_start_times and bin_size make sense
    if min([min(tt) for tt in the_list])<0: 
        raise Exception('individual_start_times values should not be negative')
    if max_value<bin_size:
        raise Exception('bin size is bigger than max individual_start_times')
    
    #specify the individual_start_times bins 
    the_bins=np.arange(0, max_value+bin_size, bin_size) # this could potentially be broken depending on what the bin_size is
    
    if the_bins[-1] < max_value:
        raise Exception('Your largest individual_start_times bin is less than max individual_start_times.  Tweak your bin size.') 
    
    middle_of_bins=np.mean(np.array([np.append(the_bins, 0), np.append(0, the_bins)]), axis=0)[1:-1]  #individual_start_times bin is defined by middle of bin
    data_in_bins=[np.array([]) for ii in range(len(middle_of_bins))] #initialize a data structure to receive data in bins
    
    #assign data to correct individual_start_times bins
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
    '''takes list of individual_start_times arrays and corresponding list of data arrays and returns the average
    by placing the data in individual_start_times bins
    '''
    data,time_bins, time_bin_middle=bin_data(time_list, data_list, bin_size)
    average_data=[]
    std_err_data=[]
    for bins in data:
        average_data.append(np.mean(bins))
        std_err_data.append(stats.sem(bins))
    assert len(average_data)==len(time_bin_middle), "data length doesn't match individual_start_times length"
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
        time_vs_psp_plot = pg.plot(labels={'left': 'individual_first_peaks of synaptic deflection (V)', 
                            'bottom': 'individual_start_times since first recorded synapse (s)', 
                            'top':(title_str+' connections: progression of synaptic defection over an experiment')})    
        ave_psp_plot = pg.plot(labels={'top':('average individual_first_baselines-line subtracted first pulse synaptic deflection ('+ title_str+ ')'), 
                          'bottom': 'individual_start_times (s)', 
                          'left':'voltage (V)'}) 
        sweep_vs_psp_plot = pg.plot(labels={'left': 'individual_first_peaks of synaptic deflection (V)', 
                            'bottom': 'sweep number', 
                            'top':(title_str+' connections: progression of synaptic defection over an experiment')})  
        
        raw=[]
        filtered=[]
        time_list=[]
        sweep_number_list=[]
        PSPs_amp_start=[]
        PSPs_amp_ave=[]
        PSPs_amp_end=[]
        length_of_experiment=[]

        slopes={}
        slopes['individual_start_times']={}
        slopes['sweep_numbers']={}
        slopes['individual_start_times']['correct_amp_good_HP']=[]
        slopes['individual_start_times']['correct_amp_bad_HP']=[]
        slopes['individual_start_times']['wrong_amp_good_HP']=[]
        slopes['individual_start_times']['wrong_amp_bad_HP']=[]
        slopes['sweep_numbers']['correct_amp_good_HP']=[]
        slopes['sweep_numbers']['correct_amp_bad_HP']=[]
        slopes['sweep_numbers']['wrong_amp_good_HP']=[]
        slopes['sweep_numbers']['wrong_amp_bad_HP']=[]
        intercepts={}
        intercepts['individual_start_times']={}
        intercepts['sweep_numbers']={}
        intercepts['individual_start_times']['correct_amp_good_HP']=[]
        intercepts['individual_start_times']['correct_amp_bad_HP']=[]
        intercepts['individual_start_times']['wrong_amp_good_HP']=[]
        intercepts['individual_start_times']['wrong_amp_bad_HP']=[]
        intercepts['sweep_numbers']['correct_amp_good_HP']=[]
        intercepts['sweep_numbers']['correct_amp_bad_HP']=[]
        intercepts['sweep_numbers']['wrong_amp_good_HP']=[]
        intercepts['sweep_numbers']['wrong_amp_bad_HP']=[]
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
    
            # figure out whether the trough or individual_first_peaks of the average synaptic trace is bigger and if that corresponds to the excitation of the neurons.  
            # i.e. if it is an excitatory synapse we would expect the max defection to be positive
            average = amp_responses.bsub_mean() #returns average synaptic response with the average baseline subtracted
            max_peak = measure_amp(average, min_or_max='max', baseline=(6e-3, 8e-3), response=(12e-3, 16e-3))
            min_peak = measure_amp(average, min_or_max='min', baseline=(6e-3, 8e-3), response=(12e-3, 16e-3))
            ave_deflection=max(abs(max_peak), abs(min_peak)) #individual_first_peaks of average trace
            max_min = "max" if abs(max_peak)> abs(min_peak) else "min"  #find whether the individual_first_peaks or trough of the first pulse average is larger 
            correct_syn_amp_dir = True
            if max_min == "min" and expt.cells[pre_id].cre_type in EXCITATORY_CRE_TYPES: 
                print ("Whoa this synapse looks inhibitory when cre line would say it should be excitatory!!!" )  
                correct_syn_amp_dir = False
            if max_min == "max" and expt.cells[pre_id].cre_type in INHIBITORY_CRE_TYPES: 
                print ("Whoa this synapse looks excitatory when cre line would say it should be inhibitory!!!" )    
                correct_syn_amp_dir = False  
       
            # find the individual_first_peaks or trough of every potential event and plot their amplitude over individual_start_times of the experiment
            individual_first_peaks=[]
            individual_first_baselines=[]
            individual_start_times=[]
            individual_holding_potentials=[]
            sweep_numbers=[]
            ordered=sorted(amp_responses.responses, key=lambda rr:rr.start_time) #order the traces by individual_start_times during the experiment
            for jj, rr in enumerate(ordered):
                individual_first_peaks.append(measure_amp(rr, min_or_max=max_min, baseline=(6e-3, 8e-3), response=(12e-3, 16e-3)))
                individual_first_baselines.append(measure_amp(rr, min_or_max=max_min, baseline=(0e-3, 2e-3), response=(6e-3, 10e-3)))
                individual_start_times.append(rr.start_time)  
                individual_holding_potentials.append(rr.parent.holding_potential)
                sweep_numbers.append(float(rr.parent.parent._sweep_id))
                print ('for each first pulse of a synapse: individual_first_peaks', individual_first_peaks[-1], 'individual_first_baselines', individual_first_baselines[-1], 'individual_start_times', individual_start_times[-1], 'holding potential', individual_holding_potentials[-1], 'sweep_numbers', sweep_numbers[-1])      

    #        for trace in amp_responses.responses:
    #            plt.plot(trace.time_values, trace.data)
    #            plt.title('responses')
    #        plt.figure()
    #        for trace in amp_responses.baselines:
    #            plt.plot(trace.time_values, trace.data)
    #            plt.title('baselines')         
    #        plt.show()
    
            # check if holding potential is within a desired range
            holding=np.mean(individual_holding_potentials) # average holding potential across plots
            print ('holding potential is', holding)
            if holding>-0.072 and holding<-0.068:
                holding_good_flag=True
            else:
                holding_good_flag=False
            print ('\tholding potential flag set to ', holding_good_flag)

            mean_base=np.mean(individual_first_baselines) # average individual_first_baselines across pulses of a synapse
            sweep_numbers=np.array(sweep_numbers)
            individual_start_times=np.array(individual_start_times) - individual_start_times[0] # remap individual_start_times basis to be in reference to start of experiment
            #TODO: !!!!from what I can tell, individual_first_peaks are already subtracting there individual baselines so I am not sure why it is being resubtracted below!!!
            peak_minus_base_average=np.array(individual_first_peaks)-mean_base # take each individual_first_peaks and put it in reference to the average individual_first_baselines
            smoothed=ndi.gaussian_filter(peak_minus_base_average, 2) # 
            t_slope, t_intercept, _,_,_=stats.linregress(individual_start_times, peak_minus_base_average)
            sn_slope, sn_intercept, _,_,_=stats.linregress(sweep_numbers, peak_minus_base_average)
            
            def update_plots_and_slopes(qc_key):
                '''updates and values for different qc groupings
                inputs: string
                    options: 'correct_amp_good_HP','correct_amp_bad_HP', 'wrong_amp_good_HP', 'wrong_amp_bad_HP'
                '''
                ave_psp_plot.plot(average.time_values, average.data, pen=pg.mkPen(color=colors[qc_key])) #plot average of first pulse in each epoch of spikes of individual synapses
                time_vs_psp_plot.plot(individual_start_times, smoothed, pen=pg.mkPen(color=colors[qc_key])) # (i, len(synapses)*1.3))
                time_vs_psp_plot.plot(individual_start_times, t_slope*individual_start_times+t_intercept, pen=pg.mkPen(color=colors[qc_key], style=pg.QtCore.Qt.DashLine)) # (i, len(synapses)*1.3))
                sweep_vs_psp_plot.plot(sweep_numbers, smoothed, pen=pg.mkPen(color=colors[qc_key]))
                sweep_vs_psp_plot.plot(sweep_numbers, sn_slope*sweep_numbers+sn_intercept, pen=pg.mkPen(color=colors[qc_key],style=pg.QtCore.Qt.DashLine)) # (i, len(synapses)*1.3))
                slopes['sweep_numbers'][qc_key].append(sn_slope)
                slopes['individual_start_times'][qc_key].append(t_slope)
                intercepts['sweep_numbers'][qc_key].append(sn_intercept)
                intercepts['individual_start_times'][qc_key].append(t_intercept)               
            # record values for different qc states
            if correct_syn_amp_dir == True and holding_good_flag ==True:
                print('recording synapse')
                time_list.append(individual_start_times)  
                filtered.append(smoothed)
                raw.append(peak_minus_base_average)
                sweep_number_list.append(sweep_numbers)
                num_of_synapses=num_of_synapses+1
                update_plots_and_slopes('correct_amp_good_HP')
                PSPs_amp_start.append(peak_minus_base_average[0])
                PSPs_amp_ave.append(np.mean(peak_minus_base_average))
                PSPs_amp_end.append(peak_minus_base_average[-1])
                length_of_experiment.append(individual_start_times[-1])
         
            else: 
                if correct_syn_amp_dir==True and holding_good_flag==False:

                    update_plots_and_slopes('correct_amp_bad_HP')
                elif correct_syn_amp_dir==False and holding_good_flag==True:
                    update_plots_and_slopes('wrong_amp_good_HP')

                elif correct_syn_amp_dir==False and holding_good_flag==False:
                    update_plots_and_slopes('wrong_amp_bad_HP')
                else:
                    print(correct_syn_amp_dir)
                    print(holding_good_flag)
                    raise Exception("This flag combo doesn't exist")
            
            print('done with one synapse')
            app.processEvents()

        
        #because times of events aren't all at the same individual_start_times, individual_start_times binning is needed to get average individual_start_times course
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
                               'sweep_std_err':sweep_std_err,
                               'slopes':slopes,
                               'intercepts':intercepts,
                               'PSPs_amp_start':PSPs_amp_start,
                               'PSPs_amp_ave':PSPs_amp_ave,
                               'PSPs_amp_end':PSPs_amp_end,
                               'length_of_experiment':length_of_experiment}
        
    ju.write("PSP_vs_time_output_data/goodpsp_vs_time or_sweep_2_07_18.json", dictionary)

    plt.figure()
    for key in dictionary.keys():
        plt.errorbar(dictionary[key]['time_points'], dictionary[key]['time_avg_data'],  yerr=dictionary[key]['time_std_err'], label=key+', n='+str(dictionary[key]['num_of_synapses']))
    plt.title('average individual_first_baselines-line subtracted first pulse synaptic deflection')
    plt.legend(loc=4)
    plt.ylabel('voltage (V)')
    plt.xlabel('individual_start_times since first recorded synapse (s)')

    plt.figure()
    for key in dictionary.keys():
        plt.errorbar(dictionary[key]['sweeps'], dictionary[key]['sweep_avg_data'],  yerr=dictionary[key]['sweep_std_err'], label=key+', n='+str(dictionary[key]['num_of_synapses']))
    plt.title('average individual_first_baselines-line subtracted first pulse synaptic deflection')
    plt.legend(loc=4)
    plt.ylabel('voltage (V)')
    plt.xlabel('sweep number')    
    
    
    plt.show(block=False)

#        app.processEvents()    
#        pg.QtGui.QApplication.exec_()  
    
    plt.show()
        
    
  
    