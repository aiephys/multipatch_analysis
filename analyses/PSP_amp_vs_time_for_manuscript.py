"""
Creates data for the elife manuscript rundown analysis using single psp
fitting but does calculate measured max for comparison.  Uses only
data pre qc'd by Steph.  Data is saved to folder specified in __main__
which should be specified for individual computers.  Data can be analyzed
using the analyze_rundown.ipynb Jupyter notebook.  This code uses 
the connection_summary not the database.
"""

from __future__ import print_function, division
import sys
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import scipy.stats as stats
import statsmodels.api as sm
from neuroanalysis.baseline import float_mode
from multipatch_analysis.connection_detection import MultiPatchSyncRecAnalyzer, MultiPatchExperimentAnalyzer
from multipatch_analysis.experiment_list import cached_experiments
from neuroanalysis.data import TraceList, PatchClampRecording
from neuroanalysis.fitting import fit_psp

# id's Steph qc'd for paper
Stephs_data=np.array([
    [('unknown', 'unknown'),	('1501090950.86', 8, 1)],
    [('unknown', 'unknown'),	('1501101571.17', 1, 5)],
    [('unknown', 'unknown'),	('1501101571.17', 1, 7)],
    [('unknown', 'unknown'),	('1501101571.17', 7, 5)],
    [('unknown', 'unknown'),	('1501104688.89', 7, 3)],
    [('unknown', 'unknown'),	('1501621744.85', 1, 6)],
    [('unknown', 'unknown'),	('1501621744.85', 6, 1)],
    [('unknown', 'unknown'),	('1501627688.56', 3, 8)],
    [('unknown', 'unknown'),	('1501627688.56', 4, 7)],
    [('unknown', 'unknown'),	('1501627688.56', 8, 3)],
    [('unknown', 'unknown'),	('1501792378.34', 2, 8)],
    [('unknown', 'unknown'),	('1501792378.34', 8, 2)],
    [('rorb', 'rorb'),	('1498687063.99', 7, 1)],
    [('rorb', 'rorb'),	('1502301827.80', 6, 8)],
    [('rorb', 'rorb'),	('1502301827.80', 8, 6)],
    [('rorb', 'rorb'),	('1523470754.85', 3, 4)],
    [('rorb', 'rorb'),	('1523470754.85', 4, 3)],
    [('rorb', 'rorb'),	('1523470754.85', 4, 6)],
    [('rorb', 'rorb'),	('1523470754.85', 4, 7)],
    [('rorb', 'rorb'),	('1523470754.85', 6, 4)], #unknown in file
    [('rorb', 'rorb'),	('1523470754.85', 7, 3)],
    [('rorb', 'rorb'),	('1523470754.85', 7, 4)],
    [('rorb', 'rorb'),	('1523470754.85', 7, 6)],
    [('rorb', 'rorb'),	('1523479910.95', 2, 3)],
    [('sim1', 'sim1'),	('1487107236.82', 7, 5)],
    [('sim1', 'sim1'),	('1487107236.82', 7, 2)],
    [('sim1', 'sim1'),	('1487367784.96', 6, 2)],
    [('sim1', 'sim1'),	('1487376645.68', 1, 7)],
    [('sim1', 'sim1'),	('1490642434.41', 5, 3)],
    [('sim1', 'sim1'),	('1490642434.41', 3, 5)],
    [('sim1', 'sim1'),	('1490642434.41', 7, 3)],
    [('sim1', 'sim1'),	('1490651407.27', 2, 5)],
    [('sim1', 'sim1'),	('1490651901.46', 4, 8)],
    [('sim1', 'sim1'),	('1497468556.18', 8, 2)],
    [('sim1', 'sim1'),	('1497468556.18', 8, 3)],
    [('sim1', 'sim1'),	('1497468556.18', 8, 6)],
    [('sim1', 'sim1'),	('1497468556.18', 2, 8)],
    [('sim1', 'sim1'),	('1497469151.70', 1, 2)],
    [('sim1', 'sim1'),	('1497469151.70', 1, 8)],
    [('sim1', 'sim1'),	('1497469151.70', 8, 5)],
    [('sim1', 'sim1'),	('1497469151.70', 8, 1)],
    [('sim1', 'sim1'),	('1497473076.69', 7, 4)],
    [('tlx3', 'tlx3'),	('1485904693.10', 8, 2)],
    [('tlx3', 'tlx3'),	('1492460382.78', 6, 2)],
    [('tlx3', 'tlx3'),	('1492460382.78', 4, 6)],
    [('tlx3', 'tlx3'),	('1492468194.97', 6, 5)],
    [('tlx3', 'tlx3'),	('1492545925.15', 2, 4)],
    [('tlx3', 'tlx3'),	('1492545925.15', 8, 5)],
    [('tlx3', 'tlx3'),	('1492545925.15', 4, 2)],
    [('tlx3', 'tlx3'),	('1492545925.15', 8, 6)],
    [('tlx3', 'tlx3'),	('1492546902.92', 2, 6)],
    [('tlx3', 'tlx3'),	('1492546902.92', 2, 8)],
    [('tlx3', 'tlx3'),	('1492546902.92', 4, 8)],
    [('tlx3', 'tlx3'),	('1492546902.92', 8, 2)],
    [('tlx3', 'tlx3'),	('1492637310.55', 5, 4)],
    [('tlx3', 'tlx3'),	('1492810479.48', 1, 7)],
    [('tlx3', 'tlx3'),	('1492812013.49', 5, 3)],
    [('tlx3', 'tlx3'),	('1494881995.55', 7, 1)],
    [('tlx3', 'tlx3'),	('1502920642.09', 7, 8)],
    [('ntsr1', 'ntsr1'),('1504737622.52', 8, 2)],
    [('ntsr1', 'ntsr1'),('1529443918.26', 1, 6)]
    ])


Steph_uids=[l[1] for l in Stephs_data]

def measure_amp_single(first_pulse_dict):
    '''measures the max of a trace within a window of the predicted psp.
    
    Input
    -----
    first_pulse_dict: dictionary that must contain following items:
        response: TraceView object
            voltage waveform of the recorded psp
        baseline: TraceView object
            voltage waveform of the baseline region
        pulse_ind: int
            index of the pulse relative to the sweep (not the local TraceView)
        rec_start: int
            index of the start of the response TraceView in reference to the original sweep data 
    
    Returns
    -------
        relative_amp: float
            measured amplitude relative to baseline
        baseline: float
            value of baseline
    '''
    response_trace=first_pulse_dict['response'].copy(t0=0) #reset time traces so can use fixed xoffset from average fit
    dt = response_trace.dt
    pulse_ind=first_pulse_dict['pulse_ind']-first_pulse_dict['rec_start'] #get the indicies of the pulse in reference to the 'respons' waveform

    psp_region_start_ind=pulse_ind+int(3e-3/dt)
    psp_region_end_ind=pulse_ind+int(15e-3/dt)

    baseline = first_pulse_dict['baseline'].mean()  #baseline voltage
    # get the maximum value in a region around where psp amp should be
    max_first_peak_region = response_trace.data[psp_region_start_ind:psp_region_end_ind].max()
    # get the minimum value in a region around where psp amp should be
    min_first_peak_region = response_trace.data[psp_region_start_ind:psp_region_end_ind].min()
    # subtract the baseline value from min and max values
    bsub=np.array([max_first_peak_region-baseline, min_first_peak_region-baseline])
    # find the absolute maximum deflection from base_line
    relative_amp=bsub[np.where(abs(bsub)==max(abs(bsub)))[0][0]] 
    # if plot==True:
    #     plt.figure()
    #     plt.baseline(response_trace.data)
    #     plt.show()
    return relative_amp, baseline


def remove_baseline_instabilities(pulse_list, baseline=2):
    """This is a qc step that Steph uses. Note: not actually used in used in final rundown analysis
    but could potentially be utilized.
    """
    for_std=[]
    for pulse in pulse_list:
        bl=pulse['baseline'].copy(data=pulse['baseline'].data-float_mode(pulse['baseline'].data), t0=0)
        for_std.append(bl)
    
    trace_mean=TraceList(for_std).mean()
    base_std=np.std(trace_mean.data)

#    base_std = np.std(for_std)
    stable_baseline=[]
    for pulse in pulse_list:
        bl=pulse['baseline'].data
        if np.abs(np.mean(bl-float_mode(bl))) < (baseline * base_std):
            stable_baseline.append(pulse)

    return stable_baseline

class fit_first_pulse():
    '''Group of functions for fitting psps of individual pulses.
    '''
    def __init__(self, expt, pre_syn_electrode_id, post_syn_electrode_id):    
        """
        Inputs
        ------
        pre_syn_electrode_id: int
            electrode id recording the pre synaptic cell
        post_syn_electrode_id: 
            electrode id recording the post synaptic cell
        pre_pad: float
            amount of time (s) before a spike to be included in the waveform.
        post_pad: float 
            amount of time (s) after a spike to be included in the waveform.
        """

        self.expt=expt
        self.pre_syn_electrode_id = pre_syn_electrode_id
        self.post_syn_electrode_id = post_syn_electrode_id
        self.pre_pad = 10e-3
        self.post_pad = 50e-3

    def get_spike_aligned_first_pulses(self):
        """Get all the first pulses that are recorded in current clamp
        and have a holding potential of between -65 and -75 and aligns 
        them at the spike times.
        
        Returns
        -------
        first_pulse_list: list of dictionaries containing information from  
            MultiPatchSyncRecAnalyzer.get.get_spike_responses() including
            the standard 'response' and 'baseline' TraceView Objects. The 
            the additional items below are added to to output dictionary.
                'sweep_id': int
                    sweep number
                'global_spike_date_time': datetime.datetime object
                    real world time of a spike. i.e the spike happend on 
                    Sept 1, 2017 at 12:34:63 pm. This is used to get the 
                    amount of time between any spike and another event such
                    as is done in 'global_seconds'
                'holding_potential': float
                    holding potential of post synaptic neuron 
                'stim_type': float
                    stimulus induction frequency in pre synaptic neuron
                'global_seconds': float
                    seconds past in reference to the first used spike in the experiment
        """
            
        # loop though sweeps in recording and pull out first pulses, in current clamp, with a holding potential between -75 mV and -65 mV
        first_pulse_list=[]
        command_trace_list=[]
        for sweep_rec in self.expt.data.contents:
            sweep_id=sweep_rec._sweep_id
            print ("SWEEP ID: %d, %s, electrode ids %d, %d, devices: %s" % (sweep_id, expt.uid, pre_syn_electrode_id, post_syn_electrode_id, sweep_rec.devices))
            if pre_syn_electrode_id not in sweep_rec.devices or post_syn_electrode_id not in sweep_rec.devices:
                print("Skipping %s electrode ids %d, %d; pre or post synaptic electrode id is not in sweep_rec.devices" % (expt.uid, pre_syn_electrode_id, post_syn_electrode_id))
                continue
            pre_rec = sweep_rec[pre_syn_electrode_id]
            post_rec = sweep_rec[post_syn_electrode_id]
            if post_rec.clamp_mode != 'ic':
                #print("Skipping %s electrode ids %d, %d; rec.clamp_mode != current clamp" % (expt.uid, pre_syn_electrode_id, post_syn_electrode_id))
                continue

            analyzer = MultiPatchSyncRecAnalyzer.get(sweep_rec)
            
            # get information about the spikes and make sure there is a spike on the first pulse
            spike_data = analyzer.get_spike_responses(pre_rec, post_rec, pre_pad=self.pre_pad, align_to='spike')
            if 0 in [pulse['pulse_n'] for pulse in spike_data]: # confirm pulse number starts at one, not zero
                raise Exception("Skipping %s electrode ids %d, %d; ; should have not have zero" % (expt.uid, pre_syn_electrode_id, post_syn_electrode_id)) 
            if 1 not in [pulse['pulse_n'] for pulse in spike_data]: # skip this sweep if not a spike on the first pulse
                #print("Skipping %s electrode ids %d, %d; no spike on first pulse" % (expt.uid, pre_syn_electrode_id, post_syn_electrode_id))                
                continue
            else: # appends sweep to the first pulse data 
                for pulse in spike_data:
                    if pulse['pulse_n'] == 1:
                        if post_rec.holding_potential<-0.075 or post_rec.holding_potential>-0.065:
                            continue     
                        pulse['sweep_id']=sweep_id
                        pulse['global_spike_date_time']=pre_rec.start_time + datetime.timedelta(0, pulse['spike']['rise_index']*pulse['response'].dt)
                        pulse['holding_potential']=post_rec.holding_potential
                        pulse['stim_type']=analyzer.stim_params(pre_rec)[0]
                        first_pulse_list.append(pulse)

        # add a key that has the time since the first used spike
        fpl_0=first_pulse_list[0]['global_spike_date_time']
        for fpl in first_pulse_list:
            fpl['global_seconds']=(fpl['global_spike_date_time']-fpl_0).total_seconds()
        
        return first_pulse_list


    def get_baseline_sub_average(self, first_pulse_list):
        """Substract the baseline for each individual fit and then
        take the average.  
        
        Input
        -----
        first_pulse_list: list of dictionaries containing information from  
            MultiPatchSyncRecAnalyzer.get.get_spike_responses() calculated
            via fit_first_pulse.get_spike_aligned_first_pulses().  Relevant 
            items for this function are
            'response': TraceView object
                post synaptic voltage waveform
            'baseline': TraceView object
                post synaptic baseline waveform
            'command': TraceView object
                current injected into the pre synaptic neuron
        """        
        bsub_trace_list=[]
        command_trace_list=[]
        for sweep in first_pulse_list:
            sweep_trace=sweep['response']
            sweep_baseline_float_mode=float_mode(sweep['baseline'].data)
            bsub_trace_list.append(sweep_trace.copy(data=sweep_trace.data-sweep_baseline_float_mode, t0=0)) #Trace object with baseline subtracted data via float mode method. Note t0 is realigned to 0 
            command_trace_list.append(sweep['command'].copy(t0=0)) #get command traces so can see the average pulse for cross talk region estimation
        
        # take average of baseline subtracted data
        self.avg_voltage=TraceList(bsub_trace_list).mean()
        self.avg_dt=self.avg_voltage.dt
        self.avg_command=TraceList(command_trace_list).mean() # pulses are slightly different in reference spike
        return self.avg_voltage, self.avg_dt, self.avg_command 

    def fit_avg(self):
        """fit the average base_line subtracted data. Updates self with the average psp fit.
        
        Returns
        -------
        self.ave_psp_fit: lmfit.model.ModelResult
            fit of the average psp waveform
        weight: numpy.ndarray
            the weight assigned to each index of the input waveform for fitting
        """
        #weighting
        weight = np.ones(len(self.avg_voltage.data))*10.  #set everything to ten initially
        weight[int((self.pre_pad-3e-3)/self.avg_dt):int(self.pre_pad/self.avg_dt)] = 0.   #area around stim artifact note that since this is spike aligned there will be some blur in where the cross talk is
        weight[int((self.pre_pad+1e-3)/self.avg_dt):int((self.pre_pad+5e-3)/self.avg_dt)] = 30.  #area around steep PSP rise 

        self.ave_psp_fit = fit_psp(self.avg_voltage, 
                        xoffset=(self.pre_pad+2e-3, self.pre_pad, self.pre_pad+5e-3), #since these are spike aligned the psp should not happen before the spike that happens at pre_pad by definition 
                        sign='any', 
                        weight=weight) 

        return self.ave_psp_fit, weight

    def fit_single(self, first_pulse_dict):
        """Fit single psp using the average fit to contrain boundries.

        Input
        -----
        first_pulse_list: list of dictionaries containing information from  
            MultiPatchSyncRecAnalyzer.get.get_spike_responses() calculated
            via fit_first_pulse.get_spike_aligned_first_pulses().  Relevant 
            items for this function are:
            'response': TraceView 
                voltage response waveform (subset of sweep TraceView)
            'pulse_ind': int
                index of the pulse relative to the sweep (not the local TraceView)
            'rec_start': int
                index of the start of the response TraceView in reference to the original sweep data 
        """
        response_trace=first_pulse_dict['response'].copy(t0=0) #reset time traces so can use fixed xoffset from average fit
#            # weight parts of the trace during fitting
        dt = response_trace.dt
        pulse_ind=first_pulse_dict['pulse_ind']-first_pulse_dict['rec_start'] #get the indicies of the pulse in reference to the 'respons' waveform
        weight = np.ones(len(response_trace.data))*10.  #set everything to ten initially
        weight[pulse_ind:pulse_ind+int(3e-3/dt)] = 0.   #area around stim artifact
        weight[pulse_ind+int(3e-3/dt):pulse_ind+int(15e-3/dt)] = 30.  #area around steep PSP rise 
        weight[pulse_ind+int(15e-3/dt):] = 0 # give decay zero weight
        
        # fit single psps while using a small jitter for xoffset, and rise_time, and fixing decay_tau
        avg_xoffset=self.ave_psp_fit.best_values['xoffset']
        xoff_min=max(avg_xoffset-.5e-3, self.pre_pad) #do not allow minimum jitter to go below the spike in the case in which xoffset of average is at the spike (which is at the location of pre_pad)
        single_psp_fit_small_bounds = fit_psp(response_trace, 
                                        xoffset=(avg_xoffset, xoff_min, avg_xoffset+.5e-3),
                                        rise_time=(self.ave_psp_fit.best_values['rise_time'], 0., self.ave_psp_fit.best_values['rise_time']),
                                        decay_tau=(self.ave_psp_fit.best_values['decay_tau'], 'fixed'),
                                        sign='any', 
                                        weight=weight)
        return single_psp_fit_small_bounds, weight



    def plot_fit(self, fit, voltage, command, weight, title, 
                measured_baseline=0,
                measured_amp=False,
                nrmse=False,
                show_plot=False,
                save_name=False):
        """plots the fit
        """

        # plot average fit
        plt.figure(figsize=(14,14))
        c1=plt.subplot(2,1,1)
        ln1=c1.plot(command.time_values*1.e3, command.data*1e3, label='command')
        c1.set_ylabel('current injection (nA)')
        c2=c1.twinx()
        ln2=c2.plot(voltage.time_values*1.e3, weight, 'k', label='weight')
        c2.set_ylabel('weight')
        plt.title('uid %s, pre/post electrodes %d, %d, individual sweeps: %s' % (self.expt.uid, 
                                    self.pre_syn_electrode_id, 
                                    self.post_syn_electrode_id, 
                                    len(first_pulse_list)))

        lines_plot_1 = ln1+ln2
        label_plot_1 = [l.get_label() for l in lines_plot_1]
        c1.legend(lines_plot_1, label_plot_1)

        ax1=plt.subplot(2,1,2)
        ln3=ax1.plot(voltage.time_values*1.e3, (voltage.data-measured_baseline)*1.e3, label='data')
        ln4=ax1.plot(voltage.time_values*1.e3, (fit.best_fit-measured_baseline)*1.e3, 'g', lw=3, label='fit')
        ax1.set_ylabel('voltage (mV)')
        ax2=ax1.twinx()
        ln5=ax2.plot(voltage.time_values*1.e3, weight, 'k', label='weight')
        ax2.set_ylabel('weight')
        ax1.set_xlabel('time (ms)')
        lines_plot_2 = ln3+ln4+ln5
        label_plot_2 = [l.get_label() for l in lines_plot_2]
        ax1.legend(lines_plot_2, label_plot_2)
        if measured_amp:
            ax1.plot(voltage.time_values*1.e3, np.ones(len(voltage.time_values))*measured_amp*1.e3, 'r--')
            plt.title(title + ' nrmse=%.3g, fit amp:%.3g, measured amp:%3g' % (nrmse, fit.best_values['amp']*1.e3, measured_amp*1.e3))
        else:
            plt.title(title + ' nrmse=%.3g, fit amp:%.3g' % (fit.nrmse(), fit.best_values['amp']*1e3))

        plt.tight_layout()
        if show_plot:
            plt.show()
        if save_name:
            plt.savefig(save_name)
            plt.close()




if __name__ == '__main__':

    # note this path must be specified for individual users
    path='/home/corinnet/workspace/aiephys/rundown_results'
    if not os.path.exists(path):
        os.makedirs(path)

    # load experiments
    expts = cached_experiments()
    dictionary={}
    synapses = []
    # cycle through all expt in connections summary (note the experiment summary is being cycled though a lot here)
    for connection in expts.connection_summary(): # expts.connection_summary() is list of dictionaries
        cells = connection['cells'] #(pre, post) synaptic "Cell" objects
        expt = connection['expt'] #"Experiment" object
        pre_syn_cell_id=cells[0].cell_id
        post_syn_cell_id=cells[1].cell_id
        pre_syn_electrode_id=cells[0].electrode.device_id
        post_syn_electrode_id=cells[1].electrode.device_id

        # skip connection if not in Stephs set 
        if (expt.uid, pre_syn_cell_id, post_syn_cell_id) not in Steph_uids:
            continue
        else:
            print ("RUNNING: %s, cell ids:%s %s, electrode ids: %s %s" % (expt.uid, pre_syn_cell_id, post_syn_cell_id, pre_syn_electrode_id, post_syn_electrode_id))

        # initialize fitting class
        fitting=fit_first_pulse(expt, pre_syn_electrode_id, post_syn_electrode_id)
        # get spike aligned first pulses at -70 holding potential
        first_pulse_list=fitting.get_spike_aligned_first_pulses()
        
        # skip neuron if there is no data
        if not len(first_pulse_list)>0:
            continue

        #this could be used to impliment Stephs curvy baseline qc.  eed to rename output above to match input here.
#        first_pulse_list=remove_baseline_instabilities(non_qc_first_pulse_list)
 #       print('original_'+ str(len(non_qc_first_pulse_list))+'_reduced_'+ str(len(first_pulse_list))+' ,'+ str(len(non_qc_first_pulse_list)-len(first_pulse_list))+"pulses were removed for baseline")


        # Get the average of the baseline subtracted first pulses
        avg_voltage, dt, avg_command=fitting.get_baseline_sub_average(first_pulse_list)
        ave_psp_fit, weight_for_average=fitting.fit_avg()

        # set up output directory and naming convention output files
        name_string=expt.uid+'_'+str(pre_syn_cell_id)+'_'+str(post_syn_cell_id)+'_'+cells[0].cre_type+'_'+cells[0].cre_type
        connection_path=os.path.join(path, name_string) #path to connection specific directory        
        if not os.path.exists(connection_path):
            os.makedirs(connection_path)

        # fit the average of the first pulses
        fitting.plot_fit(ave_psp_fit, 
                        avg_voltage, 
                        avg_command, 
                        weight_for_average, 
                        'Mean baseline subtracted spike aligned,', 
                        save_name=os.path.join(connection_path, 'AVG_'+name_string+'.png'))

        #fit individual pulses
        out_data=[]
        for first_pulse_dict in first_pulse_list:
            # fit single pulse
            single_psp_fit, weight_for_single=fitting.fit_single(first_pulse_dict)
            # get measured baseline and individual psp amplitude
            measured_amp, baseline_value=measure_amp_single(first_pulse_dict)
            
            # plot the fit and compare with measured
            fitting.plot_fit(single_psp_fit, 
                            first_pulse_dict['response'], 
                            first_pulse_dict['command'], 
                            weight_for_single, 
                            'SWEEP:'+str(first_pulse_dict['sweep_id'])+', ', 
                            nrmse=single_psp_fit.nrmse(),
                            measured_baseline=baseline_value, 
                            measured_amp=measured_amp,
                            save_name=os.path.join(connection_path, 'Sweep_'+str(first_pulse_dict['sweep_id'])+'_'+name_string+'.png'))
            
            # put data in output list
            out_data.append([expt.uid, 
                            pre_syn_cell_id, 
                            post_syn_cell_id,
                            first_pulse_dict['sweep_id'],
                            measured_amp,
                            single_psp_fit.best_values['amp'],
                            single_psp_fit.best_values['rise_time'],
                            single_psp_fit.best_values['xoffset'],
                            first_pulse_dict['global_seconds'],
                            single_psp_fit.nrmse(),
                            cells[0].cre_type,
                            cells[1].cre_type,
                            first_pulse_dict['stim_type'],
                            first_pulse_dict['global_spike_date_time']]) 
    
        # create a pandas dataframe
        out_df=pd.DataFrame(out_data)
        
        # name the columns of the dataframe
        out_df.columns=['uid',
                        'pre_cell_id',
                        'post_cell_id', 
                        'sweep_id', 
                        'measured_amp', 
                        'fit_amp', 
                        'fit_rise_time',
                        'fit_xoffset',
                        'time',
                        'nrmse', 
                        'pre_cre', 
                        'post_cre', 
                        'stim_type',
                        'global_spike_date_time']
        # save dataframe to a csv                
        out_df.to_csv(os.path.join(connection_path,'rundown.csv')) 

            
