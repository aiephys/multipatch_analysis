"""THIS IS ONLY PARTAIALLY CONVERTED"""

from neuroanalysis.data import Trace, TraceList
from multipatch_analysis.database import database as db
import multipatch_analysis.connection_strength as cs 
from multipatch_analysis.database.database import TableGroup
import matplotlib.pyplot as plt

class FirstPulseFitTableGroup(TableGroup):
    """Fits first pulse for each individual sweeps.
    """
    schemas = {

        'average_first_pulse_fit': [
            """Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed qc in current clamp""",
            ('pair_id', 'pair.id', '', {'index': True}),
            ('ic_fit_amp', 'float', 'amplitude from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_fit_latency', 'float', 'latency from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_fit_rise_time', 'float', 'rise time from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_fit_decay_tau', 'float', 'decay tau from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_amp_cv', 'float', 'coefficient of variation for first pulse amplitude in 10, 20, 50 Hz stimuli'),
            ('avg_psp', 'array', 'average psp time series, spike aligned, baseline subtracted'),
            ('n_sweeps', 'int', 'number of sweeps in avg_psp'),
            ('pulse_ids', 'object', 'list of first pulse ids in avg_psp, len(pulse_ids) == n_sweeps')
            #('ic_fit_NRMSE', 'float', 'NRMSE returned from psp_fit')
            #TODO: consider removing 50Hz responses from decay calculation
        ],
        'individual_first_pulse_fit': [
            """Best parameters fit to individual first pulses with initial conditions set 
            by the fit of the average first pulse.""",
            ('pulse_response_id', 'pulse_response.id', '', {'index': True, 'unique': True}),
#            ('pos_amp', 'float', 'max-median offset from baseline to pulse response window'),
#            ('neg_amp', 'float', 'min-median offset from baseline to pulse response window'),
#            ('pos_dec_amp', 'float', 'max-median offset from baseline to pulse response window from devonvolved trace'),
#            ('neg_dec_amp', 'float', 'min-median offset from baseline to pulse response window from deconvolved trace'),
#            ('pos_dec_latency', 'float', 'duration (seconds) from presynaptic spike max dv/dt until the sample measured in pos_dec_amp'),
#            ('neg_dec_latency', 'float', 'duration (seconds) from presynaptic spike max dv/dt until the sample measured in neg_dec_amp'),
#            ('crosstalk', 'float', 'trace difference immediately before and after onset of presynaptic stimulus pulse'),
        ],

    }
    
    def create_mappings(self):
        TableGroup.create_mappings(self)
        
        IndividualFirstPulseFits = self['individual_first_pulse_fit']
        AverageFirstPulseFits = self['average_first_pulse_fit']
        
        db.Pair.average_first_pulse_fit = db.relationship(AverageFirstPulseFits, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        AverageFirstPulseFits.pair = db.relationship(db.Pair, back_populates="average_first_pulse_fit", single_parent=True)


first_pulse_fit_tables = FirstPulseFitTableGroup()

def init_tables():
    global IndividualFirstPulseFits, AverageFirstPulseFits
    first_pulse_fit_tables.create_tables()

    IndividualFirstPulseFits = first_pulse_fit_tables['individual_first_pulse_fit']
    AverageFirstPulseFits = first_pulse_fit_tables['average_first_pulse_fit']


# create tables in database and add global variables for ORM classes
init_tables()


def update_fit(limit=None,parallel=True, workers=6, raise_exceptions=False, session=None):
    """Update table
    """
    # experiments = session.query(db.Experiment.acq_timestamp).all()
    # expts_done = session.query(db.Experiment.acq_timestamp).join(db.SyncRec).join(db.Recording).join(db.PulseResponse).join(PulseResponseStrength).distinct().all()
    # print("Skipping %d already complete experiments" % (len(expts_done)))
    # experiments = [e for e in experiments if e not in set(expts_done)]

    # get Pair objects from the database
    pairs=get_pairs_from_DB(limit)

    if parallel:
        pool = multiprocessing.Pool(processes=workers)
        pool.map(pair, pairs)
    else:
        for ii, pair in enumerate(pairs):
            # START HERE: LOOKS LIKE CONNECTION STRENGTH DOESNT HAVE STUFF
            # CHECK THAT THE WAY STEPH IS GETTING MEAN MAKES SENSE by looking at plots
            pulse_responses, pulse_ids, pulse_response_amps_measured,freq = extract_first_pulse_info_from_Pair_object(pair)
                # if len(pulse_responses) > 0:
                #     #do the fit here
                #     results = first_pulse_features(pair, pulse_responses, pulse_response_amps_measured)
                #     #format the table #TODO finish this
                #     fpf = FirstPulseFeatures(pair=pair, n_sweeps=len(pulse_ids), pulse_ids=pulse_ids, **results)
                #     s.add(fpf)
                #     if i % 10 == 0:
                #         s.commit()
                #         print("%d pairs added to the DB of %d" %(i, len(records)))
            if len(pulse_responses)>0:
                ave_psp=TraceList(pulse_responses).mean()
                plt.plot(ave_psp.data)
                plt.show()
            else:
                print ('%ii pair empty' % ii)
            print(pulse_ids)


def get_pairs_from_DB(limit):
    """Gra pairs from the database.
    input
    -----
    limit: None or integer
        specifies how many database entries to load 
    returns
    -------
    list of multipatch_analysis.database.database.Pair objects
    """
    s = db.Session() #start a database session
    #TODO: dont think need the line below (depricate when confirmed that line below works when there is something in AverageFirstPulseFits)
    #q = s.query(db.Pair, AverageFirstPulseFits).outerjoin(AverageFirstPulseFits).filter(AverageFirstPulseFits.pair_id == None)

    # get all pairs that are in pair table but are not in AverageFirstPulseFits table
    q = s.query(db.Pair, cs.ConnectionStrength).join(cs.ConnectionStrength).filter(AverageFirstPulseFits.pair_id == None)
    q = s.query("""SELECT * FROM db.Pair""")

    # perform the query with the given limits
    if limit is not None:
        q = q.limit(limit)
    print("Updating %d pairs.." % q.count())
    records=q.all()
#    return records  #return the results of the query
    # TODO: depricate line below when TODO above is resolved
    # remove irrelevant None returned from database query
    return [record[0] for record in records]


def extract_first_pulse_info_from_Pair_object(pair):
    """Extract first pulse responses and relevant information 
    from entry in the pair database. Screen out pulses that are
    not current clamp or do not pass the corresponding
    inhibitory or excitatory qc.
    input
    -----
    pair: multipatch_analysis.database.database.Pair object
    """
    # TODO: learn how to do what's below in one query
    # s = db.Session()
    # q = s.query(db.PulseResponse.data, db.StimSpike, db.PatchClampRecording)
    # q = q.join(db.StimPulse).join(db.StimSpike).join(db.PatchClampRecording)
    # filters = [
    #     (db.Pair == pair)
    #     (db.StimPulse.pulse_number == 1),
    #     (db.StimPulse.n_spikes == 1),
    #     (db.StimSpike.max_dvdt_time != None),
    #     (db.PulseResponse.ex_qc_pass == True)
    #     (db.PatchClampRecording.clamp_mode == 'ic')
    # ]
    #
    # for filter_arg in filters:
    #     q = q.filter(*filter_arg)

    try: 
        pair.connection_strength.synapse_type
    except:
        print ("skipping pair_id %s is not in connection_strength" % pair.id)
        return [], [], [], []
    synapse_type = pair.connection_strength.synapse_type
    pulse_responses = []
    pulse_response_amps_measured = []
    pulse_ids = []
    stim_freq = []
    for pr in pair.pulse_responses:
        stim_pulse = pr.stim_pulse
        n_spikes = stim_pulse.n_spikes
        pulse_number = stim_pulse.pulse_number
        pulse_id = pr.stim_pulse_id
        ex_qc_pass = pr.ex_qc_pass
        in_qc_pass = pr.in_qc_pass
        pcr = stim_pulse.recording.patch_clamp_recording
#        stim_freq = pcr.multi_patch_probe[0].induction_frequency
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
        data_trace = Trace(data=data, t0=start_time - spike_time, sample_rate=db.default_sample_rate)

        # append to output lists if neurons pass qc
        if (synapse_type == 'ex' and ex_qc_pass is True) or (synapse_type == 'in' and in_qc_pass is True):
            pulse_responses.append(data_trace)
            pulse_ids.append(pulse_id)
            stim_freq.append(pcr.multi_patch_probe[0].induction_frequency)            
        if synapse_type == 'in' and in_qc_pass is True:
            pulse_response_amps_measured.append(pr.pulse_response_strength.neg_amp)
        if synapse_type == 'ex' and ex_qc_pass is True:
            pulse_response_amps_measured.append(pr.pulse_response_strength.pos_amp)

    return pulse_responses, pulse_ids, pulse_response_amps_measured, stim_freq

if __name__=='__main__':

    update_fit(limit=None,parallel=False, workers=6, raise_exceptions=False, session=None)

