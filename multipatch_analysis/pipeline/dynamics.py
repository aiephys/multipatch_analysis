# coding: utf8
"""
For generating a DB table describing short term dynamics.

"""
from __future__ import print_function, division

import os
import numpy as np
from collections import OrderedDict
from ..util import timestamp_to_datetime
from .. import database as db
from .pipeline_module import DatabasePipelineModule
from .pulse_response import PulseResponsePipelineModule


class DynamicsPipelineModule(DatabasePipelineModule):
    """Generates dynamics analysis for each pair
    """
    name = 'dynamics'
    dependencies = [PulseResponsePipelineModule]
    table_group = db.dynamics_tables
    
    @classmethod
    def create_db_entries(cls, job_id, session):
        
        # Load experiment from DB
        expt = db.experiment_from_timestamp(job_id, session=session)
        for pair in expt.pairs.values():
            if pair.synapse is False:
                continue
            recs = session.query(db.Recording).join(db.PulseResponse).join(db.Pair).filter(db.Pair.id==pair.id).all()
            pulse_amps = {}
            for rec in recs:
                q = session.query(db.PulseResponseStrength, db.PulseResponse, db.StimPulse.pulse_number, db.MultiPatchProbe.induction_frequency)
                q = q.join(db.PulseResponse).join(db.StimPulse).join(db.PatchClampRecording, db.PatchClampRecording.recording_id==db.PulseResponse.recording_id).join(db.MultiPatchProbe)
                q = q.filter(db.PulseResponse.pair_id==pair.id).filter(db.PatchClampRecording.clamp_mode=='ic').filter(db.PulseResponse.recording_id==rec.id)
                results = q.all()
                if len(results) == 0:
                    continue
                if results[0].induction_frequency != 50:
                    continue
                sign = pair.connection_strength.synapse_type
                qc = sign+'_qc_pass'
                # make sure all 12 pulses pass qc before moving on
                qc_check = [getattr(r.Pulse_response, qc) for r in results]
                if all(qc_check) is False:
                    continue
                amp_field = 'pos_dec_amp' if sign == 'ex' else 'neg_dec_amp'
                for result in results:
                    pulse_number = result.pulse_number
                    pulse_amps.setdefault(pulse_number, [])
                    pulse_amps[pulse_number].append(getattr(result.Pulse_response_strength, amp_field))
            if any(pulse_amps):
                pulse_ratio_8_1_50Hz = np.mean(pulse_amps[8]) / np.mean(pulse_amps[1])
                pulse_ratio_2_1_50Hz = np.mean(pulse_amps[2]) / np.mean(pulse_amps[1])    
                # Write new record to DB
                dynamics = db.Dynamics(pair_id=pair.id, pulse_ratio_2_1_50Hz=pulse_ratio_2_1_50Hz, pulse_ratio_8_1_50Hz=pulse_ratio_8_1_50Hz)
                session.add(dynamics)
        
    @classmethod
    def job_records(cls, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        return session.query(db.Dynamics).filter(db.Dynamics.pair_id==db.Pair.id).filter(db.Pair.experiment_id==db.Experiment.id).filter(db.Experiment.acq_timestamp.in_(job_ids)).all()
