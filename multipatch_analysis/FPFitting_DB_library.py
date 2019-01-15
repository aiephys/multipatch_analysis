"""There are many ways the the single pulses can be fit.
These functions and objects are common between all the different 
methods.
"""

from neuroanalysis.data import Trace, TraceList
from multipatch_analysis.database import database as db
import multipatch_analysis.connection_strength as cs 
from multipatch_analysis.database.database import TableGroup
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
import matplotlib.pyplot as plt
import numpy as np
# import time
from neuroanalysis.fitting import fit_psp


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
print (len(Steph_uids))

time_before_spike = 10.e-3 #time in seconds before spike to start trace waveforms

#comon inputs to the the TableGroup
common = [('pair_id', 'pair.id', '', {'index': True}),
        ('uid', 'float','timestamp attached to the experiment for ease of viewing'),
        ('pre_cell_id', 'int', 'the integer id of the pre synaptic cell (from cell table)'),
        ('post_cell_id', 'int', 'the integer id of the post synaptic cell (from cell table)'),
        ('amp', 'float', 'amplitude '),
        ('latency', 'float', 'time elapsed since the time of presynaptic spike (max dv/dt)'),
        ('rise_time', 'float', 'rise time of psp', ),
        ('decay_tau', 'float', 'decay of psp'),
        ('avg_psp', 'array', 'array of the best fit voltage waveform starting 10 ms before pre-synaptic spike'),
        ('dt', 'float', 'time step of *avg_psp* array'),
        ('n_sweeps', 'int', 'number of sweeps used in the fit'),
        ('pulse_ids', 'object', 'data base pulse ids included in fit'),
        ('distance', 'float', 'distance between pairs'),
        ('NRMSE', 'float', 'error of fit'),
        ('synapse_sign', 'str', '"ex" or "in" (also in connection strength table but here for convenience)'),
        ('measured_baseline', 'float', 'average voltage measured between 10 and 1 ms before a spike'),
        ('measured_amp', 'float', 'amplitude within a window of 0.5 ms after spike initiation (max dv/dt) until end of array specified in the pulse_response table'),
        ('connected', 'bool', 'specifies whether human thinks their is a connection')]

def extract_first_pulse_info_from_Pair_object(pair, uid, desired_clamp='ic'):
    """Extract first pulse responses and relevant information 
    from entry in the pair database. Screen out pulses that are
    not current clamp or do not pass the corresponding
    inhibitory or excitatory qc.
    
    Input
    -----
    pair: multipatch_analysis.database.database.Pair object

    Return
    ------
    pulse_responses: TraceList of spike aligned traces where the start of each trace is 10 ms before the spike 
    pulse_ids, 
    psp_amps_measured, 
    stim_freq
    """

    
    if pair.connection_strength is None:
        print ("\t\tSKIPPING: pair_id %s, uid %s, is not yielding pair.connection_strength" % (pair.id, uid))
        return [], [], [], []
    if pair.connection_strength.synapse_type is None:
        print ("\t\tSKIPPING: pair_id %s, uid %s, is not yielding pair.connection_strength.synapse_type" % (pair.id, uid))
        return [], [], [], []
    synapse_type = pair.connection_strength.synapse_type
    pulse_responses = []
    psp_amps_measured = []
    pulse_ids = []
    stim_freqs = []
    if len(pair.pulse_responses)==0:
        print ("\t\tSKIPPING: pair_id %s, uid %s, no pulse responses in pair table" % (pair.id, uid))
        return [], [], [], []
    for pr in pair.pulse_responses:
        stim_pulse = pr.stim_pulse
        n_spikes = stim_pulse.n_spikes
        pulse_number = stim_pulse.pulse_number
        pulse_id = pr.stim_pulse_id
        ex_qc_pass = pr.ex_qc_pass
        in_qc_pass = pr.in_qc_pass
        pcr = stim_pulse.recording.patch_clamp_recording
        stim_freq = pcr.multi_patch_probe[0].induction_frequency
        clamp_mode = pcr.clamp_mode
        # current clamp
        if clamp_mode != desired_clamp:
            continue
        # ensure that there was only 1 presynaptic spike
        if n_spikes != 1:
            continue
        # we only want the first pulse of the train
        if pulse_number != 1:
            continue
        # # only include frequencies up to 50Hz
        # if stim_freq >= 100:
        #     continue

        data = pr.data
        start_time = pr.start_time
        spike_time = stim_pulse.spikes[0].max_dvdt_time        
        data_trace = Trace(data=data, t0= start_time-spike_time+time_before_spike, sample_rate=db.default_sample_rate).time_slice(start=0, stop=None) #start of the data is the spike time

        
        # append to output lists if neurons pass qc
        if (synapse_type == 'ex' and ex_qc_pass is True) or (synapse_type == 'in' and in_qc_pass is True):
            pulse_responses.append(data_trace)
            pulse_ids.append(pulse_id)
            stim_freqs.append(stim_freq)        
        if synapse_type == 'in' and in_qc_pass is True:
            psp_amps_measured.append(pr.pulse_response_strength.neg_amp)
        if synapse_type == 'ex' and ex_qc_pass is True:
            psp_amps_measured.append(pr.pulse_response_strength.pos_amp)

    return pulse_responses, pulse_ids, psp_amps_measured, stim_freq

def fit_trace(waveform, synaptic_sign, clamp_mode='ic', weight=None, latency=None, latency_jitter=None, plot_show=False, plot_save_name=False, title=''):
    """
    Input
    -----
    waveform: Trace Object
        contains data to be fit
    clamp_mode: string
        'vc' denotes voltage clamp
        'ic' denotes current clamp
    synaptic_sign: str
        'ex' or 'in' specifying excitation of synapse
    plot_show: boolean 
        show plot resulting from fit if True.
    plot_save_name: False or string
        if string is supplied then save the plot to the specified path
    title: string
        title of the resulting plot
    latency: float or None
        Amount of time that has passed in reference to the time of
        the pre-synaptic spike.  Note that this value has to be transformed in
        reference to the start of the data waveform being fit within
        this function.
    latency_jitter: None or float
        Amount of jitter to allow before and after the latency value. If
        *latency* is None, this value must be none.
    weight: numpy array
        Relative weighting of different sections of the voltage array
        for fitting.  If specified it must be the same length as the 
        waveform

    Note there is a 'time_before_spike' variable in this code that is
    not passed in and therefore must be global.  It specifies the
    the amount of data before the spike that is being considered in 
    the here and in the rest of the code.  Note that the value is positive,
    i.e. if we start looking at the trace 10 ms before the spike we use
    10e-3 not -10e-3.

    Returns
    -------
    self.ave_psp_fit: lmfit.model.ModelResult
        fit of the average psp waveform
    weight: numpy.ndarray
        the weight assigned to each index of the input waveform for fitting
    """
    #weighting
    if weight is None:
        weight = np.ones(len(waveform.data))*10.  #set everything to ten initially
        weight[int((time_before_spike-3e-3)/waveform.dt):int(time_before_spike/waveform.dt)] = 0.   #area around stim artifact note that since this is spike aligned there will be some blur in where the cross talk is
        weight[int((time_before_spike+1e-3)/waveform.dt):int((time_before_spike+5e-3)/waveform.dt)] = 30.  #area around steep PSP rise 

    #converting synaptic sign from database to convention for fitting
    if synaptic_sign == 'in':
        sign = '-'
    elif synaptic_sign == 'ex':
        sign = '+'
    elif synaptic_sign == 'any':
        sign = 'any'
    else:
        raise Exception('synaptic sign not defined')
    

    if latency is None and latency_jitter:
        raise Exception('latency_jitter cannot be specified if latency is not')
    if latency_jitter:
        if latency_jitter < .1e-3:
            raise Exception('specified latency jitter less than .0e-3 may have implications for initial conditions')
        if latency < latency_jitter:
            lower_bound = time_before_spike
        else:
            lower_bound = latency + time_before_spike - latency_jitter
        xoffset=([lower_bound, latency+time_before_spike, latency+time_before_spike+latency_jitter-.1e-3], lower_bound, latency+time_before_spike+latency_jitter)  
    elif latency:
        xoffset=(latency+time_before_spike, 'fixed')
    else:
        #since these are spike aligned the psp should not happen before the spike that happens at pre_pad by definition
        xoffset=([time_before_spike+1e-3, time_before_spike+4e-3], time_before_spike, time_before_spike+5e-3)

    fit = fit_psp(waveform, 
                    clamp_mode=clamp_mode,
                    xoffset=xoffset,
                    sign=sign, 
                    weight=weight) 

    if clamp_mode == 'ic':
        scale_factor = 1.e3
        ylabel='voltage (mV)'
    if clamp_mode == 'vc':
        scale_factor = 1.e12
        ylabel='current (pA)'        

    if plot_show is True or plot_save_name:
        plt.figure(figsize=(14,10))
        ax1=plt.subplot(1,1,1)
        ln1=ax1.plot(waveform.time_values*1.e3, waveform.data*scale_factor, 'b', label='data')
        if clamp_mode == 'ic':
            ln2=ax1.plot(waveform.time_values*1.e3, fit.best_fit*scale_factor, 'r', label='nrmse=%f \namp (mV)=%f \nlatency (ms)=%f \nrise time (ms)=%f \ndecay tau=%f' % \
                                            (fit.nrmse(), \
                                            fit.best_values['amp']*scale_factor, \
                                            (fit.best_values['xoffset']-time_before_spike)*1e3, \
                                            fit.best_values['rise_time']*1e3, \
                                            fit.best_values['decay_tau']))
        elif clamp_mode == 'vc': 
            ln2=ax1.plot(waveform.time_values*1.e3, fit.best_fit*scale_factor, 'r', label='nrmse=%f \namp (pA)=%f \nlatency (ms)=%f \nrise time (ms)=%f \ndecay tau=%f' % \
                                            (fit.nrmse(), \
                                            fit.best_values['amp']*scale_factor, \
                                            (fit.best_values['xoffset']-time_before_spike)*1e3, \
                                            fit.best_values['rise_time']*1e3, \
                                            fit.best_values['decay_tau']))
        ax2=ax1.twinx()
        ln3=ax2.plot(waveform.time_values*1.e3, weight, 'k', label='weight')
        ax1.set_ylabel(ylabel)
        ax2.set_ylabel('weight')
        ax1.set_xlabel('time (ms): spike happens at 10 ms')

        lines_plot= ln1+ln2+ln3
        label_plot = [l.get_label() for l in lines_plot]
        ax1.legend(lines_plot, label_plot)

        plt.title(title)
        if plot_show is True:
            plt.show()
        if plot_save_name:
            if plot_show is True:
                raise Exception('Cannot show and save plot')
            else: 
                plt.savefig(plot_save_name)
                plt.close()

    return fit

def measure_amp(v_array, baseline_index_window, psp_amp_index_window):
    '''measures the max of a trace within a window of the predicted psp.
    
    Input
    -----
    v_array: array-like
        voltage wave-form
    baseline_index_window: list containing two integers
        integers specify index of the beginning and end of baseline window
    psp_amp_index_window: list containing two integers
        integers specify index of the beginning and end of window to look for max/min (psp amplitude)
        
    Returns
    -------
        relative_amp: float
            measured amplitude relative to baseline
        baseline: float
            value of baseline
    '''

    # get baseline
    baseline = v_array[baseline_index_window[0]: baseline_index_window[1]].mean()  #baseline voltage
    # get the maximum value in a region around where psp amp should be
    max_psp_region = v_array[psp_amp_index_window[0]: psp_amp_index_window[1]].max()
    # get the minimum value in a region around where psp amp should be
    min_psp_region = v_array[psp_amp_index_window[0]: psp_amp_index_window[1]].min()
    # subtract the baseline value from min and max values
    bsub=np.array([max_psp_region-baseline, min_psp_region-baseline])
    # find the absolute maximum deflection from base_line
    relative_amp=bsub[np.where(abs(bsub)==max(abs(bsub)))[0][0]] 

    return relative_amp, baseline