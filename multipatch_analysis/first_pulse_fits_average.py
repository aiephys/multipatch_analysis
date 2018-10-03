#WARNING THE FORCED SIGN AND ANY SIGN SHOULD BE PULLED APART

#import pdb
from neuroanalysis.data import Trace, TraceList
from multipatch_analysis.database import database as db
import multipatch_analysis.connection_strength as cs 
from multipatch_analysis.database.database import TableGroup
import matplotlib.pyplot as plt
import numpy as np
import time
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

class FirstPulseFitTableGroup(TableGroup):
    """Fits first pulse for each individual sweeps.
    """
    schemas = {

        'avg_first_pulse_fit_any_sign': [
            """Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed qc in current clamp. Created via first_pulse_fits_average.py.
             All units in SI.""",
            ('pair_id', 'pair.id', '', {'index': True}),
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
            ('connected', 'bool', 'specifies whether human thinks their is a connection')
        ],
        'avg_first_pulse_fit_force_sign': [
            """Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed qc in current clamp. Created via first_pulse_fits_average.py.
             All units in SI.""",
            ('pair_id', 'pair.id', '', {'index': True}),
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
            ('connected', 'bool', 'specifies whether human thinks their is a connection')
        ]
        # 'individual_first_pulse_fit': [
        #     """Best parameters fit to individual first pulses with initial conditions set 
        #     by the fit of the average first pulse available in the avg_first_pulse_fit_force_sign table.""",
        #     ('pulse_response_id', 'pulse_response.id', '', {'index': True, 'unique': True}),
        #     ('pair_id', 'pair.id', '', {'index': True}),
        #     ('amp', 'float', ''),
        #     ('latency', 'float', '(seconds) from presynaptic spike (max dv/dt)'),
        #     ('rise_time', 'float', ''),
        #     ('decay_tau', 'float', ''),
        #     ('avg_psp', 'array', ''),
        #     ('n_sweeps', 'int', ''),
        #     ('pulse_ids', 'object', ''),
        #     ('NRMSE', 'float', '')
        # ]

    }
    
    def create_mappings(self):
        TableGroup.create_mappings(self)
        
#        IndividualFirstPulseFits = self['individual_first_pulse_fit']
        AvgFirstPulseFitsAnySign = self['avg_first_pulse_fit_any_sign']
        AvgFirstPulseFitsForceSign = self['avg_first_pulse_fit_force_sign']
        
        db.Pair.avg_first_pulse_fit_force_sign = db.relationship(AvgFirstPulseFitsForceSign, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        AvgFirstPulseFitsForceSign.pair = db.relationship(db.Pair, back_populates="avg_first_pulse_fit_force_sign", single_parent=True)
        
        db.Pair.avg_first_pulse_fit_any_sign = db.relationship(AvgFirstPulseFitsAnySign, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        AvgFirstPulseFitsAnySign.pair = db.relationship(db.Pair, back_populates="avg_first_pulse_fit_any_sign", single_parent=True)


first_pulse_fit_tables = FirstPulseFitTableGroup()

def init_tables():
    global AvgFirstPulseFitsAnySign
    global AvgFirstPulseFitsForceSign 
#    global IndividualFirstPulseFits, 
    first_pulse_fit_tables.create_tables()

#    IndividualFirstPulseFits = first_pulse_fit_tables['individual_first_pulse_fit']
    AvgFirstPulseFitsAnySign = first_pulse_fit_tables['avg_first_pulse_fit_any_sign']
    AvgFirstPulseFitsForceSign = first_pulse_fit_tables['avg_first_pulse_fit_force_sign']


# create tables in database and add global variables for ORM classes
init_tables()

def update_fit(limit=None, expts=None, parallel=True, workers=6, raise_exceptions=False, session=None):
    """
    """
    session=db.Session()
    if expts is None:
        experiments = session.query(db.Experiment.acq_timestamp).all()
        #TODO: confirm this query is good enough
        expts_done=session.query(db.Experiment.acq_timestamp).join(db.Pair).join(AvgFirstPulseFitsAnySign).all()
#        #TODO: deprecate this line when line above is confirmed expts_done = session.query(db.Experiment.acq_timestamp).join(db.SyncRec).join(db.Recording).join(AverageFirstPulseFits).join().distinct().all()
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

def compute_fit(job_info, raise_exceptions=False):
    
    session = db.Session() #create session

    expt_id, index, n_jobs = job_info
    print("QUERYING (expt_id=%f): %d/%d" % (expt_id, index, n_jobs))

    #do the query
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

    processed_count = 0 #index for keeping track of how many cells pairs in experiemnt have been analized
    for ii, (pair, uid, pre_cell_id, post_cell_id, pre_cell_cre, post_cell_cre) in enumerate(expt_stuff):
    #            skip connection if not in Stephs set 
    #            if (str(np.round(uid, 2)), pre_cell_id, post_cell_id) not in Steph_uids:
    #                #print ("SKIPPING: %s, cell ids:%s %s" % (uid, pre_cell_id, post_cell_id))
    #                continue
    #            else:
        print ("\tTRYING TO GET FIRST PULSES: number %i of %i experiment pairs: %s, cell ids:%s %s" % (ii, len(expt_stuff), uid, pre_cell_id, post_cell_id))
        pulse_responses, pulse_ids, psp_amps_measured, freq = extract_first_pulse_info_from_Pair_object(pair, uid)
        
        # if pulses are returned fit the average
        if len(pulse_responses)>0:
            print ("\tFITTING: %s, cell ids:%s %s" % (uid, pre_cell_id, post_cell_id))

            avg_psp=TraceList(pulse_responses).mean()
    #                for pr in pulse_responses:
    #                    plt.plot(pr.time_values, pr.data)
    #                plt.plot(ave_psp.time_values, ave_psp.data, lw=5)
    #                plt.show()
        else:
            print ("\t\tSKIPPING: %s, cell ids:%s %s: no passing pulse responses" % (uid, pre_cell_id, post_cell_id))                                                           
            continue

        # get the measured baseline and amplitude of psp
        measured_relative_amp, measured_baseline=measure_amp(avg_psp.data, 
                            [0, int((time_before_spike-1.e-3)/avg_psp.dt)], 
                            [int((time_before_spike+.5e-3)/avg_psp.dt), -1])

        # deal with when there is not distance measurement in pair table
        if pair.distance:
            pair_distance=pair.distance
        else: pair_distance =np.infty 

        # grab syapse from the table
        synapse_sign=pair.connection_strength.synapse_type
        connected=pair.synapse

        # specify where to save files
        title='%s, cells %i %s to %i %s; distance=%.1f um' % (uid, pre_cell_id, pre_cell_cre, post_cell_id,post_cell_cre, pair_distance*1e6)
        save_name='/home/corinnet/workspace/DBfit_pics_any_sign/%s_%s%s_%s%s_average_fit.png'  % (uid, pre_cell_id, pre_cell_cre, post_cell_id,post_cell_cre)
        
        # fit the trace
        avg_fit=fit_trace(avg_psp, synaptic_sign='any', plot_save_name=save_name, title=title)
        
        # specify results and fit them
        result_dict={'dt' : avg_psp.dt,
                    'amp': avg_fit.best_values['amp'], 
                    'latency': avg_fit.best_values['xoffset']-time_before_spike,
                    'rise_time':  avg_fit.best_values['rise_time'],
                    'decay_tau': avg_fit.best_values['decay_tau'],
                    'avg_psp': avg_fit.best_fit,
                    'NRMSE': avg_fit.nrmse()}

        afpf=AvgFirstPulseFitsAnySign(pair=pair, 
                                    distance=pair_distance, 
                                    uid=uid, 
                                    pre_cell_id=pre_cell_id,
                                    post_cell_id=post_cell_id,
                                    n_sweeps=len(pulse_ids), 
                                    pulse_ids=pulse_ids,
                                    synapse_sign=synapse_sign, 
                                    measured_amp=measured_relative_amp,
                                    measured_baseline=measured_baseline,
                                    connected=connected,
                                    **result_dict)
        session.add(afpf)
        processed_count=processed_count+1

            
        #     # s.add(afpf)
        #     # if i % 10 == 0:
        #     #     s.commit()
        #     #     print("%d pairs added to the DB of %d" %(i, len(records)))              
    session.commit()
    print("COMMITED %i pairs from expt_id=%f: %d/%d" % (processed_count, expt_id, index, n_jobs))

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

def fit_trace(voltage_trace, synaptic_sign, plot_show=False, plot_save_name=False, title=''):
    """
    Input
    -----
    voltage_trace: Trace Object
        contains data to be fit
    synaptic_sign: str
        'ex' or 'in' specifying excitation of synapse
    plot_show: boolean 
        plot resulting fit if True
    plot_save_name: False or string
        if string is supplied then save the plot to the specified path
    title: string
        title of the resulting plot

    Returns
    -------
    self.ave_psp_fit: lmfit.model.ModelResult
        fit of the average psp waveform
    weight: numpy.ndarray
        the weight assigned to each index of the input waveform for fitting
    """
    #weighting
    weight = np.ones(len(voltage_trace.data))*10.  #set everything to ten initially
    weight[int((time_before_spike-3e-3)/voltage_trace.dt):int(time_before_spike/voltage_trace.dt)] = 0.   #area around stim artifact note that since this is spike aligned there will be some blur in where the cross talk is
    weight[int((time_before_spike+1e-3)/voltage_trace.dt):int((time_before_spike+5e-3)/voltage_trace.dt)] = 30.  #area around steep PSP rise 

    #converting synaptic sign from database to convention for fitting
    if synaptic_sign == 'in':
        sign = '-'
    elif synaptic_sign == 'ex':
        sign = '+'
    elif synaptic_sign == 'any':
        sign = 'any'
    else:
        raise Exception('synaptic sign not defined')


    fit = fit_psp(voltage_trace, 
                    xoffset=(time_before_spike+2e-3, time_before_spike, time_before_spike+5e-3), #since these are spike aligned the psp should not happen before the spike that happens at pre_pad by definition 
                    sign=sign, 
                    weight=weight) 

    if plot_show is True or plot_save_name:
        plt.figure(figsize=(14,14))
        ax1=plt.subplot(1,1,1)
        ln1=ax1.plot(voltage_trace.time_values*1.e3, voltage_trace.data, 'b', label='data')
        ln2=ax1.plot(voltage_trace.time_values*1.e3, fit.best_fit, 'r', label='nrmse=%f' % fit.nrmse())
        ax2=ax1.twinx()
        ln3=ax2.plot(voltage_trace.time_values*1.e3,weight, 'k', label='weight')
        ax1.set_ylabel('voltage')
        ax2.set_ylabel('weight')

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

def extract_first_pulse_info_from_Pair_object(pair, uid):
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
        if clamp_mode != 'ic':
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




if __name__=='__main__':

    #Note that after this is done being prototyped delete so dont accedently overwrite table
#    first_pulse_fit_tables.drop_tables()
    init_tables()
#    update_fit(limit=None, expts=[1533768797.736], parallel=False, workers=6, raise_exceptions=False, session=None)
    update_fit(limit=None, expts=[1492545925.146], parallel=False, workers=6, raise_exceptions=False, session=None)

#    update_fit(limit=None, expts=None, parallel=False, workers=6, raise_exceptions=False, session=None)
