"""Fit average of first pulses
"""

from neuroanalysis.data import Trace, TraceList
from multipatch_analysis.database import database as db
import multipatch_analysis.connection_strength as cs 
from multipatch_analysis.database.database import TableGroup
import matplotlib.pyplot as plt
import numpy as np
import time
from neuroanalysis.fitting import fit_psp
#import FPFitting_DB_library as FPF_lib
import datetime
import os

commiting = True
time_before_spike = 10.e-3 #time in seconds before spike to start trace waveforms

class FirstPulseFitTableGroup(TableGroup):
    """Fits first pulse for each individual sweeps.
    """
    schemas = {
        'avg_first_pulse_fit': [
             """Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed qc in current clamp. During the fit the latency is forced
            to be the value +/-.5 ms found via fitting all of the pulses (available in the 
            connection_strength.ic_fit_xoffset). The heavily weighted section (meant to 
            place more importance of the wave form during the rise time) is shifted to 
            begin at the latency. Created via fit_average_first_pulse.py. 
            All units in SI.""",
            ('pair_id', 'pair.id', 'The ID of the entry in the pair table to which these results apply', {'index': True}),

            # current clamp
            ('ic_amp', 'float', 'fit amplitude of current clamp average first pulses'),
            ('ic_latency', 'float', 'fit time elapsed since the time of presynaptic spike (max dv/dt) of current clamp data'),
            ('ic_rise_time', 'float', 'fit rise time of psp of current clamp data'),
            ('ic_decay_tau', 'float', 'fit decay of psp of current clamp data'),
            ('ic_avg_psp', 'array', 'fit array of the best fit voltage waveform starting 10 ms before pre-synaptic spike'),
            ('ic_dt', 'float', 'time step of *avg_psp* array from current clamp data'),
            ('ic_pulse_ids', 'object', 'data base pulse ids included in the current clamp fit'),
            ('ic_NRMSE', 'float', 'error of fit of current clamp fit'),
            ('ic_measured_baseline', 'float', 'average voltage measured between 10 and 1 ms before a spike'),
            ('ic_measured_amp', 'float', 'voltage amplitude within a window of 0.5 ms after spike initiation (max dv/dt) until end of array specified in the pulse_response table'),
            ('ic_weight', 'array', 'weighting used during fitting of current clamp data'),

            # voltage clamp
            ('vc_amp', 'float', 'fit amplitude of voltage clamp average first pulses'),
            ('vc_latency', 'float', 'fit time elapsed since the time of presynaptic spike (max dv/dt) of voltage clamp data'),
            ('vc_rise_time', 'float', 'fit rise time of psp measured in voltage clamp'),
            ('vc_decay_tau', 'float', 'fit decay of psp measured in voltage clamp'),
            ('vc_avg_psp', 'array', 'fit array of the best fit current waveform starting 10 ms before pre-synaptic spike'),
            ('vc_dt', 'float', 'time step of *avg_psp* array from voltage clamp data'),
            ('vc_pulse_ids', 'object', 'data base pulse ids included in the voltage clamp fit'),
            ('vc_NRMSE', 'float', 'error of fit of voltage clamp fit'),
            ('vc_measured_baseline', 'float', 'average current measured between 10 and 1 ms before a spike'),
            ('vc_measured_amp', 'float', 'current amplitude within a window of 0.5 ms after spike initiation (max dv/dt) until end of array specified in the pulse_response table'),
            ('vc_weight', 'array', 'weighting used during fitting of voltage clamp data')]
    }

    def create_mappings(self):
        TableGroup.create_mappings(self)

        AvgFirstPulseFit = self['avg_first_pulse_fit']
        db.Pair.avg_first_pulse_fit = db.relationship(AvgFirstPulseFit, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        AvgFirstPulseFit.pair = db.relationship(db.Pair, back_populates="avg_first_pulse_fit", single_parent=True)

first_pulse_fit_tables = FirstPulseFitTableGroup()

def init_tables():

    global AvgFirstPulseFit
    first_pulse_fit_tables.create_tables()
    AvgFirstPulseFit = first_pulse_fit_tables['avg_first_pulse_fit']

# create tables in database and add global variables for ORM classes
init_tables()

def update_DB(limit=None, expts=None, parallel=True, workers=6, raise_exceptions=False, session=None):
    """
    """
    session=db.Session()
    if expts is None:
        experiments = session.query(db.Experiment.acq_timestamp).all()
        expts_done=session.query(db.Experiment.acq_timestamp).join(db.Pair).join(AvgFirstPulseFit).all()
        print("Skipping %d already complete experiments" % (len(expts_done)))
        experiments = [e for e in experiments if e not in set(expts_done)]

        if limit > 0:
            np.random.shuffle(experiments)
            experiments = experiments[:limit]

        jobs = [(record.acq_timestamp, index, len(experiments)) for index, record in enumerate(experiments)]
    else:
        jobs = [(expt, i, len(expts)) for i, expt in enumerate(expts)]
    # if parallel:
    #     pool = multiprocessing.Pool(processes=workers)
    #     pool.map(pair, pairs)
    # else:
    for job in jobs:
        compute_fit(job, raise_exceptions=raise_exceptions)

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

def extract_first_pulse_info_from_Pair_object(pair, desired_clamp='ic'):
    """Extract first pulse responses and relevant information 
    from entry in the pair database. Screen out pulses that are
    not current clamp or do not pass the corresponding
    inhibitory or excitatory qc.
    
    Input
    -----
    pair: multipatch_analysis.database.database.Pair object
    desired_clamp: string
        Specifies whether current or voltage clamp sweeps are desired.
        Options are:
            'ic': current clamp
            'vc': voltage clamp
    
    Return
    ------
    pulse_responses: TraceList 
        traces where the start of each trace is 10 ms before the spike 
    pulse_ids: list of ints
        pulse ids of *pulse_responses*
    psp_amps_measured: list of floats
        amplitude of *pulse_responses* from the *pulse_response* table
    stim_freq: list of floats
        the stimulation frequency corresponding to the *pulse_responses* 
    """

    if pair.connection_strength is None:
        print ("\t\tSKIPPING: pair_id %s, is not yielding pair.connection_strength" % pair.id)
        return [], [], [], []
    if pair.connection_strength.synapse_type is None:
        print ("\t\tSKIPPING: pair_id %s, is not yielding pair.connection_strength.synapse_type" % pair.id)
        return [], [], [], []
    synapse_type = pair.connection_strength.synapse_type
    pulse_responses = []
    psp_amps_measured = []
    pulse_ids = []
    stim_freqs = []
    if len(pair.pulse_responses)==0:
        print ("\t\tSKIPPING: pair_id, no pulse responses in pair table" % (pair.id))
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


def get_average_pulse_response(pair, desired_clamp='ic'):
    """
    Inputs
    ------
    pair: multipatch_analysis.database.database.Pair object

    desired_clamp: string
        Specifies whether current or voltage clamp sweeps are desired.
        Options are:
            'ic': current clamp
            'vc': voltage clamp

    Returns
    -------
    pulse_responses
    pulse_ids
    psp_amps_measured
    freq
    avg_psp
    measured_relative_amp
    measured_baseline
    """
    # get pulses that pass qc
    pulse_responses, pulse_ids, psp_amps_measured, freq = extract_first_pulse_info_from_Pair_object(pair, desired_clamp=desired_clamp)

    # if pulses are returned take the average
    if len(pulse_responses)>0:
        avg_psp=TraceList(pulse_responses).mean()
#                for pr in pulse_responses:
#                    plt.plot(pr.time_values, pr.data)
#                plt.plot(ave_psp.time_values, ave_psp.data, lw=5)
#                plt.show()
    else:
        return None, None, None, None, None, None, None

    # get the measured baseline and amplitude of psp
    measured_relative_amp, measured_baseline=measure_amp(avg_psp.data, 
                        [0, int((time_before_spike-1.e-3)/avg_psp.dt)], 
                        [int((time_before_spike+.5e-3)/avg_psp.dt), -1])

    return pulse_responses, pulse_ids, psp_amps_measured, freq, avg_psp, measured_relative_amp, measured_baseline 

def fit_trace(waveform, excitation, clamp_mode='ic', weight=None, latency=None, latency_jitter=None):
    """
    Input
    -----
    waveform: Trace Object
        contains data to be fit
    clamp_mode: string
        'vc' denotes voltage clamp
        'ic' denotes current clamp
    excitation: str
        'ex' or 'in' specifying excitation of synapse
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
    # weighting
    if weight is None:
        weight = np.ones(len(waveform.data))  #set everything to ones

    # set fitting sign to positive or negative based on excitation and clamp state
    if (excitation == 'in') and (clamp_mode == 'ic'):
        sign = '-'
    elif (excitation == 'in') and (clamp_mode == 'vc'):
        sign = '+'
    elif (excitation == 'ex') and (clamp_mode == 'ic'):
        sign = '+'
    elif (excitation == 'ex') and (clamp_mode == 'vc'):
        sign = '-'
    elif excitation == 'any':
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

    if False:
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
        plt.show()

    return fit


def compute_fit(job_info, raise_exceptions=False):
    
    session = db.Session() #create session

    expt_id, index, n_jobs = job_info
    print("QUERYING (expt_id=%f): %d/%d" % (expt_id, index, n_jobs))

    #do query
    pre_cell = db.aliased(db.Cell)
    post_cell = db.aliased(db.Cell)
    expt_stuff = session.query(db.Pair, db.Experiment.acq_timestamp, pre_cell.ext_id, post_cell.ext_id,pre_cell.cre_type, post_cell.cre_type)\
                        .join(db.Experiment)\
                        .join(pre_cell, db.Pair.pre_cell_id==pre_cell.id)\
                        .join(post_cell, db.Pair.post_cell_id==post_cell.id).filter(db.Experiment.acq_timestamp==expt_id).all()
    # make sure query returned something
    if len(expt_stuff) <=0:
        print('No pairs found for expt_id=%f', expt_id)
        return

    processed_count = 0 #index for keeping track of how many cells pairs in experiemnt have been analyzed
    for ii, (pair, uid, pre_cell_id, post_cell_id, pre_cell_cre, post_cell_cre) in enumerate(expt_stuff):

        print ("Number %i of %i experiment pairs: %0.3f, cell ids:%s %s" % (ii, len(expt_stuff), uid, pre_cell_id, post_cell_id))
        
        # grab syapse from the table
        try:
            excitation=pair.connection_strength.synapse_type
        except:
            print('\tskipping: no pair.connection_strength.synapse_type')
            continue

        if not pair.connection_strength.ic_fit_xoffset:
            print('\tskipping: no latency to do forced latency fitting')
            continue
        xoffset=pair.connection_strength.ic_fit_xoffset

        # -----------fit current clamp data---------------------        
        # get pulses
        (pulse_responses_i, pulse_ids_i, psp_amps_measured_i, freq, avg_psp_i, 
            measured_relative_amp_i, measured_baseline_i) = get_average_pulse_response(pair, desired_clamp='ic')

        if pulse_responses_i:
            # weight and fit the trace
            weight_i = np.ones(len(avg_psp_i.data))*10.  #set everything to ten initially
            weight_i[int((time_before_spike-3e-3)/avg_psp_i.dt):int(time_before_spike/avg_psp_i.dt)] = 0.   #area around stim artifact note that since this is spike aligned there will be some blur in where the cross talk is
            weight_i[int((time_before_spike+.0001+xoffset)/avg_psp_i.dt):int((time_before_spike+.0001+xoffset+4e-3)/avg_psp_i.dt)] = 30.  #area around steep PSP rise 
            avg_fit_i = fit_trace(avg_psp_i, excitation=excitation, weight=weight_i, latency=xoffset, latency_jitter=.5e-3)
            latency_i = avg_fit_i.best_values['xoffset']-time_before_spike
            amp_i = avg_fit_i.best_values['amp']
            rise_time_i = avg_fit_i.best_values['rise_time']
            decay_tau_i = avg_fit_i.best_values['decay_tau']
            avg_fit_waveform_i = avg_fit_i.best_fit
            dt_i = avg_psp_i.dt
            nrmse_i = avg_fit_i.nrmse()
        else:
            print('\tskipping: no suitable first pulses found in current clamp')
            weight_i = None
            latency_i = None
            amp_i = None
            rise_time_i = None
            decay_tau_i = None
            avg_fit_waveform_i = None
            dt_i = None
            nrmse_i = None
        # --------------fit voltage clamp data---------------------        
        # get pulses
        (pulse_responses_v, pulse_ids_v, psp_amps_measured_v, freq_v, avg_psp_v,  
            measured_relative_amp_v, measured_baseline_v) = get_average_pulse_response(pair, desired_clamp='vc')

        if pulse_responses_v:
            # weight and fit the trace    
            weight_v = np.ones(len(avg_psp_v.data))*10.  #set everything to ten initially
            weight_v[int((time_before_spike+.0001+xoffset)/avg_psp_v.dt):int((time_before_spike+.0001+xoffset+4e-3)/avg_psp_v.dt)] = 30.  #area around steep PSP rise 
            avg_fit_v = fit_trace(avg_psp_v, excitation=excitation, clamp_mode = 'vc', weight=weight_v, latency=xoffset, latency_jitter=.5e-3)
            latency_v = avg_fit_v.best_values['xoffset']-time_before_spike
            amp_v = avg_fit_v.best_values['amp']
            rise_time_v = avg_fit_v.best_values['rise_time']
            decay_tau_v = avg_fit_v.best_values['decay_tau']
            avg_fit_waveform_v = avg_fit_v.best_fit
            dt_v = avg_psp_v.dt
            nrmse_v = avg_fit_v.nrmse()

        else:
            print('\tskipping: no suitable first pulses found in voltage clamp')
            weight_v = None
            latency_v = None
            amp_v = None
            rise_time_v = None
            decay_tau_v = None
            avg_fit_waveform_v = None
            dt_v = None
            nrmse_v = None
        #------------ done with fitting section ------------------------------

        # dictionary for ease of translation into the output table
        out_dict={
             'ic_amp': amp_i,
             'ic_latency': latency_i,
             'ic_rise_time': rise_time_i,
             'ic_decay_tau': decay_tau_i,
             'ic_avg_psp': avg_fit_waveform_i,
             'ic_dt': dt_i,
             'ic_pulse_ids': pulse_ids_i,
             'ic_NRMSE': nrmse_i,
             'ic_measured_baseline': measured_baseline_i,
             'ic_measured_amp': measured_relative_amp_i,
             'ic_weight': np.array(weight_i),

             'vc_amp': amp_v,
             'vc_latency': latency_v,
             'vc_rise_time': rise_time_v,
             'vc_decay_tau': decay_tau_v,
             'vc_avg_psp': avg_fit_waveform_v,
             'vc_dt': dt_v,
             'vc_pulse_ids': pulse_ids_v,
             'vc_NRMSE': nrmse_v,
             'vc_measured_baseline': measured_baseline_v,
             'vc_measured_amp': measured_relative_amp_v,
             'vc_weight': np.array(weight_v)
             } 
        # map to pair table and commit
        afpf=AvgFirstPulseFit(pair=pair, **out_dict)
        if commiting is True:
            session.add(afpf)
            session.commit()
#---------------------------------------------------------------------------------------        
        processed_count=processed_count+1
        print('processed', processed_count+1)

    # if commiting is True:
    #     # pair.meta = pair.meta.copy()  # required by sqlalchemy to flag as modified
    #     # pair.meta['avg_first_pulse_fit_timestamp'] = time.time()  

        print("COMMITED %i pairs from expt_id=%f: %d/%d" % (processed_count, expt_id, index, n_jobs))

if __name__=='__main__':

#    first_pulse_fit_tables.drop_tables() #note this will drop all the tables here!
    init_tables()
#    update_DB(limit=None, expts=[1533768797.736], parallel=False, workers=6, raise_exceptions=False, session=None)

    update_DB(limit=None, expts=None, parallel=False, workers=6, raise_exceptions=False, session=None)