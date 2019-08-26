# coding: utf8
from __future__ import print_function, division

import os
from sqlalchemy.orm import aliased
from ... import config
from ..pipeline_module import DatabasePipelineModule
from .synapse_prediction import SynapsePredictionPipelineModule
from ...fit_average_first_pulse import fit_average_first_pulses, fit_single_first_pulse
import traceback
import sys


class AverageFirstPulseFitPipelineModule(DatabasePipelineModule):
    """Measure the "resting state" amplitude of a synapse by averaging together only responses 
    that have no prior stimuli in a certain window. 
    """
    name = 'avg_first_pulse_fit'
    dependencies = [SynapsePredictionPipelineModule]
    table_group = ['avg_first_pulse_fit']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        expt_id = job['job_id']

        expt = db.experiment_from_timestamp(expt_id, session=session)
       
        fails = []
        errors = ""
        for (pre_cell_id, post_cell_id), pair in expt.pairs.items():
            try:
                result = fit_average_first_pulses(pair)
                if result['error'] is not None:
                    # known error occurred; we consider this a successful run
                    errors += "(%s->%s) %s\n\n" % (pre_cell_id, post_cell_id, result['error'])
                else:
                    result.pop('error')
                    afpf = db.AvgFirstPulseFit(pair=pair, **result)
                    session.add(afpf)
            except:
                # unhandled exception occurred; we consider this an unsuccessful run
                exc = traceback.format_exception(*sys.exc_info())
                errors += "(%s->%s) Exception\n%s" % (pre_cell_id, post_cell_id, exc)
                fails.append((pre_cell_id, post_cell_id))
               
        if len(fails) > 0:
            raise Exception("Exceptions occurred processing pairs: %s\n\n%s" % (fails, errors))
        else:
            return errors
        
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        q = session.query(db.AvgFirstPulseFit)
        q = q.filter(db.AvgFirstPulseFit.pair_id==db.Pair.id)
        q = q.filter(db.Pair.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.ext_id.in_(job_ids))
        return q.all()


def join_pulse_response_to_expt(db, query):
    pre_rec = aliased(db.Recording)
    post_rec = aliased(db.Recording)
    joins = [
        (post_rec, db.PulseResponse.recording),
        (db.PatchClampRecording,),
        (db.MultiPatchProbe,),
        (db.SyncRec,),
        (db.Experiment,),
        (db.StimPulse, db.PulseResponse.stim_pulse),
        (pre_rec, db.StimPulse.recording),
    ]
    for join_args in joins:
        query = query.join(*join_args)

    return query, pre_rec, post_rec

def get_single_pulse_data(session, pair, db):
    """Select records from pulse_response_strength table
    This is a version grabbed from cs.get_amps altered to get
    additional info.  Note that if cs.get_baseline_amps is also
    called there could be some unexpected results/correspondense
    problems if cs.get_amps or cs.get_baseline is altered.
    """
    cols = [
        db.PulseResponse.id,
        db.PulseResponse.ex_qc_pass,
        db.PulseResponse.in_qc_pass,
        db.PatchClampRecording.clamp_mode,
        db.PatchClampRecording.baseline_potential,
        db.PatchClampRecording.baseline_current,
        db.StimPulse.pulse_number,
        db.StimSpike.max_slope_time.label('spike_time'), #note that in my fit_average_first_pulse code I filter for spikes[0] so there may be an error with multiple here
        db.PulseResponse.data_start_time.label('response_start_time'),
        db.MultiPatchProbe.induction_frequency.label('stim_freq'),
        db.PulseResponse.data,
    ]


    q = session.query(*cols)
    #q = q.join(db.PulseResponse)
    
    q, pre_rec, post_rec = join_pulse_response_to_expt(db, q)
    q = q.join(db.StimSpike)
    # note that these are aliased so needed to be joined outside the join_pulse_response_to_expt() function
    q = q.add_columns(post_rec.start_time.label('rec_start_time'))
    q = q.add_columns(post_rec.id.label('post_rec_id'))  #this will identify the train the pulse is in
        
    filters = [
        (pre_rec.electrode==pair.pre_cell.electrode,),
        (post_rec.electrode==pair.post_cell.electrode,),
        (db.PatchClampRecording.qc_pass==True,),
        (db.StimPulse.pulse_number==1,),
    ]
    for filter_args in filters:
        q = q.filter(*filter_args)
    
    # should result in chronological order
    q = q.order_by(db.PulseResponse.id).all()

    return(q)
    # df = pandas.read_sql_query(q.statement, q.session.bind)
    # recs = df.to_records()
    # return recs




class SingleFirstPulseFitPipelineModule(DatabasePipelineModule):
    """Analyze synaptic connection strength for all pairs per experiment
    """
    name = 'single_first_pulse_fit'
    dependencies = [SynapsePredictionPipelineModule]
    table_group = ['single_first_pulse_fit']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        expt_id = job['job_id']
        expt = db.experiment_from_timestamp(expt_id, session=session)
       
        fails = []
        errors = ""
        for (pre_cell_id, post_cell_id), pair in expt.pairs.items():
            passing_first_pulses = get_single_pulse_data(session, pair, db)         
            for pr in passing_first_pulses: 
                try:
                    result = fit_single_first_pulse(pr, pair)
                    if result['error'] is not None:
                        # known error occurred; we consider this a successful run
                        errors += "(%s->%s) %s\n\n" % (pre_cell_id, post_cell_id, result['error'])
                    else:
                        result.pop('error')
                        afpf = db.SingleFirstPulseFit(pulse_response_id=pr.id, **result)
                        session.add(afpf)
                except:
                    # unhandled exception occurred; we consider this an unsuccessful run
                    exc = traceback.format_exception(*sys.exc_info())
                    errors += "(%s->%s) Exception\n%s" % (pre_cell_id, post_cell_id, exc)
                    fails.append((pre_cell_id, post_cell_id))
               
        if len(fails) > 0:
            raise Exception("Exceptions occurred processing pairs: %s\n\n%s" % (fails, errors))
        else:
            return errors
        
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        q = session.query(db.SingleFirstPulseFit)
        q = q.filter(db.SingleFirstPulseFit.pulse_response_id==db.PulseResponse.id)
        q = q.filter(db.PulseResponse.pair_id==db.Pair.id)
        q = q.filter(db.Pair.experiment_id==db.Experiment.id)
        q = q.filter(db.Experiment.acq_timestamp.in_(job_ids))

        return q.all()