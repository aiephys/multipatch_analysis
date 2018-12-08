#WARNING THE FORCED SIGN AND ANY SIGN SHOULD BE PULLED APART
#ADD LATENCY TO THE MIX

#import pdb
from neuroanalysis.data import Trace, TraceList
from multipatch_analysis.database import database as db
import multipatch_analysis.connection_strength as cs 
from multipatch_analysis.database.database import TableGroup
import matplotlib.pyplot as plt
import numpy as np
import time
from neuroanalysis.fitting import fit_psp
import FPFitting_DB_library as FPF_lib
import datetime
import os
#----------------------------------------------------------
#----- specify things up here------------------------------
#----------------------------------------------------------

fitting_type = 'vclamp' #dynamic_weight_jitter_latency' #options 'default', 'force_sign', force_latency, force_sign_and_latency, dynamic_weight, dynamic_weight_jitter_latency
general_image_path='/home/corinnet/workspace/DBfit_pics/'
save_image = True    #specifies whether to save images
commiting = True
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------

date=datetime.datetime.today().strftime('%Y-%m-%d')

if fitting_type == 'default':
    image_path=os.path.join(general_image_path, 'default'+date)
elif fitting_type == 'force_sign':
    image_path=os.path.join(general_image_path, 'force_sign'+date)
elif fitting_type == 'force_latency':
    image_path=os.path.join(general_image_path, 'force_latency'+date)
elif fitting_type == 'force_sign_and_latency':
    image_path=os.path.join(general_image_path, 'force_SL'+date)
elif fitting_type == 'dynamic_weight':
    image_path=os.path.join(general_image_path, 'dynamic_weight'+date)
elif fitting_type == 'dynamic_weight_jitter_latency':
    image_path=os.path.join(general_image_path, 'dynamic_weight_jitter_latency'+date)
elif fitting_type == 'vclamp':
    image_path=os.path.join(general_image_path, 'vclamp'+date)
else:
    raise Exception('A recognized type of fitting has not been specified')

if save_image==True:
    confirm_save = raw_input("You are preforming %s fitting. Save fitting images to %s? " % (fitting_type, image_path)) == 'y'
    if not os.path.exists(image_path):
        os.makedirs(image_path)
else:
    print('WARNING YOU ARE NOT SAVING FIGURES')

class FirstPulseFitTableGroup(TableGroup):
    """Fits first pulse for each individual sweeps.
    """
    schemas = {
        'avg_first_pulse_fit_default': [
            """Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed qc in current clamp. No constraints on fitting other than 
            latency is within the reasonable boundries of between 0 and 5 ms after the max dv/dt
            of the spike.  Created via first_pulse_fits_average.py.
             All units in SI."""]+FPF_lib.common,
        'avg_first_pulse_fit_force_sign': [
            """Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed qc in current clamp. During the fit the sign is forced
            to be in the direction specified in the connection strength table. Created via 
            first_pulse_fits_average.py. All units in SI."""]+FPF_lib.common,
        'avg_first_pulse_fit_force_latency': [
            """Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed qc in current clamp. During the fit the latency is forced
            to be the value found via fitting all of the pulses (available in  
            connection_strengthic_fit_xoffset). The sign is forced to be in the direction 
            specified in the connection strength tableCreated via first_pulse_fits_average.py. 
            All units in SI."""]+FPF_lib.common,
        'avg_first_pulse_fit_force_SL': [
             """Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed qc in current clamp. During the fit the latency is forced
            to be the value found via fitting all of the pulses (available in the 
            connection_strengthic_fit_xoffset). Created via first_pulse_fits_average.py. 
            All units in SI."""]+FPF_lib.common,
        'avg_first_pulse_fit_dynamic_weight': [
             """Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed qc in current clamp. During the fit the latency is forced
            to be the value found via fitting all of the pulses (available in the 
            connection_strength.ic_fit_xoffset). The heavily weighted section meant to 
            place more importance of the wave form during the rise time is shifted to 
            begin at the latency. Created via first_pulse_fits_average.py. 
            All units in SI."""]+FPF_lib.common,
        'avg_first_pulse_fit_dynamic_w_latency_jitter2': [
             """Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed qc in current clamp. During the fit the latency is forced
            to be the value +/-.5 ms found via fitting all of the pulses (available in the 
            connection_strength.ic_fit_xoffset). The heavily weighted section meant to 
            place more importance of the wave form during the rise time is shifted to 
            begin at the latency. Created via first_pulse_fits_average.py. 
            All units in SI."""]+FPF_lib.common,
        'avg_first_pulse_fit_vclamp4': [
            """Fits of average first pulse voltage clamp traces."""]+FPF_lib.common
    }
    
    def create_mappings(self):
        TableGroup.create_mappings(self)

#THIS IS A REPEAT EXCEPT THE TABLES NAME SPECIFICATION 
#-----------------------------------------------------------------------------------        
        AvgFirstPulseFitsDefault = self['avg_first_pulse_fit_default']
        db.Pair.avg_first_pulse_fit_default = db.relationship(AvgFirstPulseFitsDefault, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        AvgFirstPulseFitsDefault.pair = db.relationship(db.Pair, back_populates="avg_first_pulse_fit_default", single_parent=True)
#----------------------------------------------------------------------------------        
        AvgFirstPulseFitsForceSign = self['avg_first_pulse_fit_force_sign']
        db.Pair.avg_first_pulse_fit_force_sign = db.relationship(AvgFirstPulseFitsForceSign, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        AvgFirstPulseFitsForceSign.pair = db.relationship(db.Pair, back_populates="avg_first_pulse_fit_force_sign", single_parent=True)
#-----------------------------------------------------------------------------------        
        AvgFirstPulseFitsForceLatency = self['avg_first_pulse_fit_force_latency']
        db.Pair.avg_first_pulse_fit_force_latency = db.relationship(AvgFirstPulseFitsForceLatency, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        AvgFirstPulseFitsForceLatency.pair = db.relationship(db.Pair, back_populates="avg_first_pulse_fit_force_latency", single_parent=True)
#-----------------------------------------------------------------------------------
        AvgFirstPulseFitsForceSL = self['avg_first_pulse_fit_force_SL']
        db.Pair.avg_first_pulse_fit_force_SL = db.relationship(AvgFirstPulseFitsForceSL, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        AvgFirstPulseFitsForceSL.pair = db.relationship(db.Pair, back_populates="avg_first_pulse_fit_force_SL", single_parent=True)
#-----------------------------------------------------------------------------------
        AvgFirstPulseFitsDynamicWeight = self['avg_first_pulse_fit_dynamic_weight']
        db.Pair.avg_first_pulse_fit_dynamic_weight = db.relationship(AvgFirstPulseFitsDynamicWeight, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        AvgFirstPulseFitsDynamicWeight.pair = db.relationship(db.Pair, back_populates="avg_first_pulse_fit_dynamic_weight", single_parent=True)
#-----------------------------------------------------------------------------------
        AFPFForceSignJitterLatency = self['avg_first_pulse_fit_dynamic_w_latency_jitter2']
        db.Pair.avg_first_pulse_fit_dynamic_w_latency_jitter2 = db.relationship(AFPFForceSignJitterLatency, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        AFPFForceSignJitterLatency.pair = db.relationship(db.Pair, back_populates="avg_first_pulse_fit_dynamic_w_latency_jitter2", single_parent=True)
#-----------------------------------------------------------------------------------
        VClamp = self['avg_first_pulse_fit_vclamp4']
        db.Pair.avg_first_pulse_fit_vclamp4 = db.relationship(VClamp, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        VClamp.pair = db.relationship(db.Pair, back_populates="avg_first_pulse_fit_vclamp4", single_parent=True)

first_pulse_fit_tables = FirstPulseFitTableGroup()

def init_tables():
    global AvgFirstPulseFitsDefault
    global AvgFirstPulseFitsForceSign
    global AvgFirstPulseFitsForceLatency     
    global AvgFirstPulseFitsForceSL
    global AvgFirstPulseFitsDynamicWeight
    global AFPFForceSignJitterLatency
    global VClamp
#    global IndividualFirstPulseFits, 
    first_pulse_fit_tables.create_tables()

#    IndividualFirstPulseFits = first_pulse_fit_tables['individual_first_pulse_fit']
    AvgFirstPulseFitsDefault = first_pulse_fit_tables['avg_first_pulse_fit_default']
    AvgFirstPulseFitsForceSign = first_pulse_fit_tables['avg_first_pulse_fit_force_sign']
    AvgFirstPulseFitsForceLatency = first_pulse_fit_tables['avg_first_pulse_fit_force_latency']
    AvgFirstPulseFitsForceSL = first_pulse_fit_tables['avg_first_pulse_fit_force_SL']
    AvgFirstPulseFitsDynamicWeight = first_pulse_fit_tables['avg_first_pulse_fit_dynamic_weight']
    AFPFForceSignJitterLatency = first_pulse_fit_tables['avg_first_pulse_fit_dynamic_w_latency_jitter2']
    VClamp = first_pulse_fit_tables['avg_first_pulse_fit_vclamp4']

# create tables in database and add global variables for ORM classes
init_tables()

def update_fit(limit=None, expts=None, parallel=True, workers=6, raise_exceptions=False, session=None):
    """
    """
    session=db.Session()
    if expts is None:
        experiments = session.query(db.Experiment.acq_timestamp).all()
        #TODO: confirm this query is good enough
        if fitting_type == 'default':
            expts_done=session.query(db.Experiment.acq_timestamp).join(db.Pair).join(AvgFirstPulseFitsDefault).all()
        elif fitting_type == 'force_sign':
            expts_done=session.query(db.Experiment.acq_timestamp).join(db.Pair).join(AvgFirstPulseFitsForceSign).all()
        elif fitting_type == 'force_latency':
            expts_done=session.query(db.Experiment.acq_timestamp).join(db.Pair).join(AvgFirstPulseFitsForceLatency).all()
        elif fitting_type == 'force_sign_and_latency':
            expts_done=session.query(db.Experiment.acq_timestamp).join(db.Pair).join(AvgFirstPulseFitsForceSL).all()
        elif fitting_type == 'dynamic_weight':
            # this moves the heavily weighted section meant to happen during the rise time to a pre
            # defined latency found by fitting all of the psps together
            expts_done=session.query(db.Experiment.acq_timestamp).join(db.Pair).join(AvgFirstPulseFitsDynamicWeight).all() 
        elif fitting_type == 'dynamic_weight_jitter_latency':
            expts_done=session.query(db.Experiment.acq_timestamp).join(db.Pair).join(AFPFForceSignJitterLatency).all()
        elif fitting_type == 'vclamp':
            expts_done=session.query(db.Experiment.acq_timestamp).join(db.Pair).join(VClamp).all()


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
        print ("\tTRYING TO GET FIRST PULSES: number %i of %i experiment pairs: %0.3f, cell ids:%s %s" % (ii, len(expt_stuff), uid, pre_cell_id, post_cell_id))
        if fitting_type == 'vclamp':
            pulse_responses, pulse_ids, psp_amps_measured, freq = FPF_lib.extract_first_pulse_info_from_Pair_object(pair, uid, desired_clamp='vc')
        else:
            pulse_responses, pulse_ids, psp_amps_measured, freq = FPF_lib.extract_first_pulse_info_from_Pair_object(pair, uid)
        
        # if pulses are returned fit the average
        if len(pulse_responses)>0:
            print ("\tFITTING: %0.3f, cell ids:%s %s" % (uid, pre_cell_id, post_cell_id))

            avg_psp=TraceList(pulse_responses).mean()
    #                for pr in pulse_responses:
    #                    plt.plot(pr.time_values, pr.data)
    #                plt.plot(ave_psp.time_values, ave_psp.data, lw=5)
    #                plt.show()
        else:
            print ("\t\tSKIPPING: %s, cell ids:%s %s: no passing pulse responses" % (uid, pre_cell_id, post_cell_id))                                                           
            continue

        # get the measured baseline and amplitude of psp
        measured_relative_amp, measured_baseline=FPF_lib.measure_amp(avg_psp.data, 
                            [0, int((FPF_lib.time_before_spike-1.e-3)/avg_psp.dt)], 
                            [int((FPF_lib.time_before_spike+.5e-3)/avg_psp.dt), -1])

        # deal with when there is not distance measurement in pair table
        if pair.distance:
            pair_distance=pair.distance
        else: pair_distance =np.infty 

        # grab syapse from the table
        synapse_sign=pair.connection_strength.synapse_type
        connected=pair.synapse

        title='%0.3f, cells %i %s to %i %s; distance=%.1f um; n=%i' % (uid, pre_cell_id, pre_cell_cre, post_cell_id,post_cell_cre, pair_distance*1e6, len(pulse_ids))
        #Specify name for plot saving.
        if save_image:
            save_image_name=os.path.join(image_path,'%0.3f_%s%s_%s%s_average_fit.png'  % (uid, pre_cell_id, pre_cell_cre, post_cell_id,post_cell_cre))
        else:
            save_image_name=None

        # fit the trace----------------------------------------------------------------------------------
        if fitting_type == 'default':
            avg_fit = FPF_lib.fit_trace(avg_psp, synaptic_sign='any', plot_save_name=save_image_name, title=title)
        elif fitting_type == 'force_sign':
            avg_fit = FPF_lib.fit_trace(avg_psp, synaptic_sign=synapse_sign, plot_save_name=save_image_name, title=title)
        elif fitting_type == 'force_latency':
            if pair.connection_strength.ic_fit_xoffset:
                xoffset=pair.connection_strength.ic_fit_xoffset
                avg_fit = FPF_lib.fit_trace(avg_psp, synaptic_sign='any', latency=xoffset, plot_save_name=save_image_name, title=title)
            else:
                print('No latency to do forced latency fitting')
                continue
        elif fitting_type == 'force_sign_and_latency':
            if pair.connection_strength.ic_fit_xoffset:
                xoffset=pair.connection_strength.ic_fit_xoffset
                avg_fit = FPF_lib.fit_trace(avg_psp, synaptic_sign=synapse_sign, latency=xoffset, plot_save_name=save_image_name, title=title)
            else:
                print('No latency to do forced latency fitting')
                continue
        elif fitting_type == 'dynamic_weight':
            if pair.connection_strength.ic_fit_xoffset:
                xoffset=pair.connection_strength.ic_fit_xoffset
                weight = np.ones(len(avg_psp.data))*10.  #set everything to ten initially
                weight[int((FPF_lib.time_before_spike-3e-3)/avg_psp.dt):int(FPF_lib.time_before_spike/avg_psp.dt)] = 0.   #area around stim artifact note that since this is spike aligned there will be some blur in where the cross talk is
                weight[int((FPF_lib.time_before_spike+.0001+xoffset)/avg_psp.dt):int((FPF_lib.time_before_spike+.0001+xoffset+4e-3)/avg_psp.dt)] = 30.  #area around steep PSP rise 
                avg_fit = FPF_lib.fit_trace(avg_psp, synaptic_sign=synapse_sign, weight=weight, latency=xoffset, plot_save_name=save_image_name, title=title)
            else:
                print('No latency to do forced latency fitting')
                continue
        elif fitting_type == 'dynamic_weight_jitter_latency':
            if pair.connection_strength.ic_fit_xoffset:
                xoffset=pair.connection_strength.ic_fit_xoffset
                weight = np.ones(len(avg_psp.data))*10.  #set everything to ten initially
                weight[int((FPF_lib.time_before_spike-3e-3)/avg_psp.dt):int(FPF_lib.time_before_spike/avg_psp.dt)] = 0.   #area around stim artifact note that since this is spike aligned there will be some blur in where the cross talk is
                weight[int((FPF_lib.time_before_spike+.0001+xoffset)/avg_psp.dt):int((FPF_lib.time_before_spike+.0001+xoffset+4e-3)/avg_psp.dt)] = 30.  #area around steep PSP rise 
                avg_fit = FPF_lib.fit_trace(avg_psp, synaptic_sign=synapse_sign, weight=weight, latency=xoffset, latency_jitter=.5e-3, plot_save_name=save_image_name, title=title)
            else:
                print('No latency to do forced latency fitting')
                continue
        elif fitting_type == 'vclamp':
            if pair.connection_strength.ic_fit_xoffset:
                xoffset=pair.connection_strength.ic_fit_xoffset
                weight = np.ones(len(avg_psp.data))*10.  #set everything to ten initially
#                weight[int((FPF_lib.time_before_spike-3e-3)/avg_psp.dt):int(FPF_lib.time_before_spike/avg_psp.dt)] = 0.   #area around stim artifact note that since this is spike aligned there will be some blur in where the cross talk is
                weight[int((FPF_lib.time_before_spike+.0001+xoffset)/avg_psp.dt):int((FPF_lib.time_before_spike+.0001+xoffset+4e-3)/avg_psp.dt)] = 30.  #area around steep PSP rise 
                #TODO: FIX THIS IN THE CODE fit_trace FUNCTION
                if synapse_sign == 'ex':
                    synapse_sign = 'in'
                else:
                    synapse_sign = 'ex'
                avg_fit = FPF_lib.fit_trace(avg_psp, synaptic_sign=synapse_sign, clamp_mode = 'vc', weight=weight, latency=xoffset, latency_jitter=.5e-3, plot_save_name=save_image_name, title=title)
            else:
                print('No latency to do forced latency fitting')
                continue
        else:
            raise Exception('The type of fitting has not been specified')
#---------------------------------------------------------------------------------------------

        # dictionary for ease of translation into the output table
        out_dict={'dt' : avg_psp.dt,
                    'amp': avg_fit.best_values['amp'], 
                    'latency': avg_fit.best_values['xoffset']-FPF_lib.time_before_spike,
                    'rise_time':  avg_fit.best_values['rise_time'],
                    'decay_tau': avg_fit.best_values['decay_tau'],
                    'avg_psp': avg_fit.best_fit,
                    'NRMSE': avg_fit.nrmse(),
                    'distance': pair_distance, 
                    'uid': uid, 
                    'pre_cell_id': pre_cell_id,
                    'post_cell_id': post_cell_id,
                    'n_sweeps': len(pulse_ids), 
                    'pulse_ids': pulse_ids,
                    'synapse_sign': synapse_sign, 
                    'measured_amp': measured_relative_amp,
                    'measured_baseline': measured_baseline,
                    'connected': connected}

#---------------------------------------------------------------------------------------
        if fitting_type == 'default':
            afpf=AvgFirstPulseFitsDefault(pair=pair, **out_dict)
        elif fitting_type == 'force_sign':
            afpf=AvgFirstPulseFitsForceSign(pair=pair, **out_dict)
        elif fitting_type == 'force_latency':
            afpf=AvgFirstPulseFitsForceLatency(pair=pair, **out_dict)
        elif fitting_type == 'force_sign_and_latency':
            afpf=AvgFirstPulseFitsForceSL(pair=pair, **out_dict)
        elif fitting_type == 'dynamic_weight':
            afpf=AvgFirstPulseFitsDynamicWeight(pair=pair, **out_dict)
        elif fitting_type == 'dynamic_weight_jitter_latency':
            afpf=AFPFForceSignJitterLatency(pair=pair, **out_dict)
        elif fitting_type == 'vclamp':
            afpf=VClamp(pair=pair, **out_dict)
        else:
            raise Exception('The type of fitting has not been specified')
        if commiting is True:
            session.add(afpf)
#---------------------------------------------------------------------------------------        
        processed_count=processed_count+1

    if commiting is True:
        # pair.meta = pair.meta.copy()  # required by sqlalchemy to flag as modified
        # pair.meta['Corinne_timestamp'] = time.time()  
        session.commit()
        print("COMMITED %i pairs from expt_id=%f: %d/%d" % (processed_count, expt_id, index, n_jobs))

if __name__=='__main__':

    #Note that after this is done being prototyped delete so dont accedently overwrite table
#    first_pulse_fit_tables.drop_tables() #note this will drop all the tables here!
    init_tables()
#    update_fit(limit=None, expts=[1533768797.736], parallel=False, workers=6, raise_exceptions=False, session=None)
#    update_fit(limit=None, expts=[1492545925.146], parallel=False, workers=6, raise_exceptions=False, session=None)
#    update_fit(limit=None, expts=[1533765069.190], parallel=False, workers=6, raise_exceptions=False, session=None)
#    update_fit(limit=None, expts=[1496435165.254], parallel=False, workers=6, raise_exceptions=False, session=None)
#    update_fit(limit=None, expts=[1507235159.776], parallel=False, workers=6, raise_exceptions=False, session=None)
#    update_fit(limit=None, expts=[1486511976.478], parallel=False, workers=6, raise_exceptions=False, session=None)
#    update_fit(limit=None, expts=[1524091806.551], parallel=False, workers=6, raise_exceptions=False, session=None)
#    update_fit(limit=None, expts=[1486769398.428], parallel=False, workers=6, raise_exceptions=False, session=None)
#    update_fit(limit=None, expts=[1523658475.023], parallel=False, workers=6, raise_exceptions=False, session=None)
#    update_fit(limit=None, expts=[1486753952.16], parallel=False, workers=6, raise_exceptions=False, session=None)
#    update_fit(limit=None, expts=[1483648804.078], parallel=False, workers=6, raise_exceptions=False, session=None)
#    update_fit(limit=None, expts=[1499890136.257], parallel=False, workers=6, raise_exceptions=False, session=None)

    update_fit(limit=None, expts=None, parallel=False, workers=6, raise_exceptions=False, session=None)
