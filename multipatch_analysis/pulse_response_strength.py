# coding: utf8
"""
Analyses that measure the strength of synaptic connections.

"""
from __future__ import print_function, division

import sys, multiprocessing

import numpy as np
import pyqtgraph as pg

from neuroanalysis.data import Trace
from neuroanalysis import filter
from neuroanalysis.event_detection import exp_deconvolve
from neuroanalysis.baseline import float_mode

from .database import database as db
from .database import TableGroup


class PulseResponseStrengthTableGroup(TableGroup):
    """Measures pulse amplitudes for each pulse response and background chunk.
    """
    schemas = {
        'pulse_response_strength': [
            "Measurements of membrane potential or current deflection following each evoked presynaptic spike.",
            ('pulse_response_id', 'pulse_response.id', '', {'index': True, 'unique': True}),
            ('pos_amp', 'float', 'max-median offset from baseline to pulse response window'),
            ('neg_amp', 'float', 'min-median offset from baseline to pulse response window'),
            ('pos_dec_amp', 'float', 'max-median offset from baseline to pulse response window from devonvolved trace'),
            ('neg_dec_amp', 'float', 'min-median offset from baseline to pulse response window from deconvolved trace'),
            ('pos_dec_latency', 'float', 'duration (seconds) from presynaptic spike max dv/dt until the sample measured in pos_dec_amp'),
            ('neg_dec_latency', 'float', 'duration (seconds) from presynaptic spike max dv/dt until the sample measured in neg_dec_amp'),
            ('crosstalk', 'float', 'trace difference immediately before and after onset of presynaptic stimulus pulse'),
        ],
        'baseline_response_strength' : [
            "Measurements of membrane potential or current deflection in the absence of presynaptic spikes (provides a measurement of background noise to compare to pulse_response_strength).",
            ('baseline_id', 'baseline.id', '', {'index': True, 'unique': True}),
            ('pos_amp', 'float', 'max-median offset from baseline to pulse response window'),
            ('neg_amp', 'float', 'min-median offset from baseline to pulse response window'),
            ('pos_dec_amp', 'float', 'max-median offset from baseline to pulse response window from devonvolved trace'),
            ('neg_dec_amp', 'float', 'min-median offset from baseline to pulse response window from deconvolved trace'),
            ('pos_dec_latency', 'float', 'duration (seconds) from presynaptic spike max dv/dt until the sample measured in pos_dec_amp'),
            ('neg_dec_latency', 'float', 'duration (seconds) from presynaptic spike max dv/dt until the sample measured in neg_dec_amp'),
            ('crosstalk', 'float', 'trace difference immediately before and after onset of presynaptic stimulus pulse'),
        ]
    }

    def create_mappings(self):
        TableGroup.create_mappings(self)
        
        PulseResponseStrength = self['pulse_response_strength']
        BaselineResponseStrength = self['baseline_response_strength']
        
        db.PulseResponse.pulse_response_strength = db.relationship(PulseResponseStrength, back_populates="pulse_response", cascade="delete", single_parent=True, uselist=False)
        PulseResponseStrength.pulse_response = db.relationship(db.PulseResponse, back_populates="pulse_response_strength", single_parent=True)

        db.Baseline.baseline_response_strength = db.relationship(BaselineResponseStrength, back_populates="baseline", cascade="delete", single_parent=True, uselist=False)
        BaselineResponseStrength.baseline = db.relationship(db.Baseline, back_populates="baseline_response_strength", single_parent=True)




pulse_response_strength_tables = PulseResponseStrengthTableGroup()


def init_tables():
    global PulseResponseStrength, BaselineResponseStrength
    pulse_response_strength_tables.create_tables()

    PulseResponseStrength = pulse_response_strength_tables['pulse_response_strength']
    BaselineResponseStrength = pulse_response_strength_tables['baseline_response_strength']


# create tables in database and add global variables for ORM classes
init_tables()


def measure_peak(trace, sign, spike_time, pulse_times, spike_delay=1e-3, response_window=4e-3):
    # Start measuring response after the pulse has finished, and no earlier than 1 ms after spike onset
    # response_start = max(spike_time + spike_delay, pulse_times[1])

    # Start measuring after spike and hope that the pulse offset doesn't get in the way
    # (if we wait for the pulse to end, then we miss too many fast rise / short latency events)
    response_start = spike_time + spike_delay
    response_stop = response_start + response_window

    # measure baseline from beginning of data until 50µs before pulse onset
    baseline_start = 0
    baseline_stop = pulse_times[0] - 50e-6

    baseline = float_mode(trace.time_slice(baseline_start, baseline_stop).data)
    response = trace.time_slice(response_start, response_stop)

    if sign == '+':
        i = np.argmax(response.data)
    else:
        i = np.argmin(response.data)
    peak = response.data[i]
    latency = response.time_values[i] - spike_time
    return peak - baseline, latency


def measure_sum(trace, sign, baseline=(0e-3, 9e-3), response=(12e-3, 17e-3)):
    baseline = trace.time_slice(*baseline).data.sum()
    peak = trace.time_slice(*response).data.sum()
    return peak - baseline

        
def deconv_filter(trace, pulse_times, tau=15e-3, lowpass=1000., lpf=True, remove_artifacts=False, bsub=True):
    dec = exp_deconvolve(trace, tau)

    if remove_artifacts:
        # after deconvolution, the pulse causes two sharp artifacts; these
        # must be removed before LPF
        cleaned = remove_crosstalk_artifacts(dec, pulse_times)
    else:
        cleaned = dec

    if bsub:
        baseline = np.median(cleaned.time_slice(cleaned.t0+5e-3, cleaned.t0+10e-3).data)
        b_subbed = cleaned-baseline
    else:
        b_subbed = cleaned

    if lpf:
        return filter.bessel_filter(b_subbed, lowpass)
    else:
        return b_subbed


def remove_crosstalk_artifacts(data, pulse_times):
    dt = data.dt
    r = [-50e-6, 250e-6]
    edges = [(int((t+r[0])/dt), int((t+r[1])/dt)) for t in pulse_times]
    # If window is too shortm then it becomes seneitive to sample noise.
    # If window is too long, then it becomes sensitive to slower signals (like the AP following pulse onset)
    return filter.remove_artifacts(data, edges, window=100e-6)


@db.default_session
def update_strength(limit=0, expts=None, parallel=True, workers=6, raise_exceptions=False, session=None):
    """Update pulse response strength tables for all experiments
    """
    if expts is None:
        experiments = session.query(db.Experiment.acq_timestamp).all()
        expts_done = session.query(db.Experiment.acq_timestamp).join(db.SyncRec).join(db.Recording).join(db.Baseline).join(BaselineResponseStrength).distinct().all()
        print("Skipping %d already complete experiments" % (len(expts_done)))
        experiments = [e for e in experiments if e not in set(expts_done)]

        if limit > 0:
            np.random.shuffle(experiments)
            experiments = experiments[:limit]

        jobs = [(record.acq_timestamp, index, len(experiments)) for index, record in enumerate(experiments)]
    else:
        jobs = [(expt, i, len(expts)) for i, expt in enumerate(expts)]

    if parallel:
        pool = multiprocessing.Pool(processes=workers)
        pool.map(compute_strength, jobs)
    else:
        for job in jobs:
            compute_strength(job, raise_exceptions=raise_exceptions)


def compute_strength(job_info, raise_exceptions=False):
    """Fill pulse_response_strength and baseline_response_strength tables for all pulse responses in the given experiment.
    """
    session = db.Session()
    
    try:
        expt_id, index, n_jobs = job_info
        print("Analyzing pulse response strength (expt_id=%f): %d/%d" % (expt_id, index, n_jobs))
        _compute_strength('pulse_response', expt_id, session=session)
        _compute_strength('baseline', expt_id, session=session)
        session.commit()
    except:
        session.rollback()
        print("Error in experiment: %f" % expt_id)
        if raise_exceptions:
            raise
        else:
            sys.excepthook(*sys.exc_info())
    

@db.default_session
def _compute_strength(source, expt_id, session=None):
    """Compute per-pulse-response strength metrics
    """
    if source == 'baseline':
        q = baseline_query(session)
    elif source == 'pulse_response':
        q = response_query(session)
    else:
        raise ValueError("Invalid source %s" % source)

    # select just data for the selected experiment
    q = q.join(db.SyncRec).join(db.Experiment).filter(db.Experiment.acq_timestamp==expt_id)

    prof = pg.debug.Profiler(delayed=False)
    
    recs = q.all()
    prof('fetch')
        
    new_recs = []

    for rec in recs:
        result = analyze_response_strength(rec, source)
        new_rec = {'%s_id'%source: rec.response_id}
        # copy a subset of results over to new record
        for k in ['pos_amp', 'neg_amp', 'pos_dec_amp', 'neg_dec_amp', 'pos_dec_latency', 'neg_dec_latency', 'crosstalk']:
            new_rec[k] = result[k]
        new_recs.append(new_rec)
    
    prof('process')

    # Bulk insert is not safe with parallel processes
    # if source == 'pulse_response':
    #     session.bulk_insert_mappings(PulseResponseStrength, new_recs)
    # else:
    #     session.bulk_insert_mappings(BaselineResponseStrength, new_recs)
    if source == 'pulse_response':
        for rec in new_recs:
            prs = PulseResponseStrength(**rec)
            session.add(prs)
    else:
        for rec in new_recs:
            brs = BaselineResponseStrength(**rec)
            session.add(brs)

    prof('insert')
    new_recs = []

    return "succeeded"


def response_query(session):
    """
    Build a query to get all pulse responses along with presynaptic pulse and spike timing
    """
    q = session.query(
        db.PulseResponse.id.label('response_id'),
        db.PulseResponse.data,
        db.PulseResponse.start_time.label('rec_start'),
        db.StimPulse.onset_time.label('pulse_start'),
        db.StimPulse.duration.label('pulse_dur'),
        db.StimSpike.max_dvdt_time.label('spike_time'),
        db.PatchClampRecording.clamp_mode,
        db.PulseResponse.ex_qc_pass,
        db.PulseResponse.in_qc_pass,
    )
    q = q.join(db.StimPulse)
    q = q.join(db.StimSpike)
    q = q.join(db.PulseResponse.recording).join(db.PatchClampRecording) 

    # return qc-failed records as well so we can verify qc is working
    #q = q.filter(((db.PulseResponse.ex_qc_pass==True) | (db.PulseResponse.in_qc_pass==True)))

    return q


def baseline_query(session):
    """
    Build a query to get all baseline responses
    """
    q = session.query(
        db.Baseline.id.label('response_id'),
        db.Baseline.data,
        db.PatchClampRecording.clamp_mode,
        db.Baseline.ex_qc_pass,
        db.Baseline.in_qc_pass,
    )
    q = q.join(db.Baseline.recording).join(db.PatchClampRecording) 

    # return qc-failed records as well so we can verify qc is working
    # q = q.filter(((db.Baseline.ex_qc_pass==True) | (db.Baseline.in_qc_pass==True)))

    return q


def analyze_response_strength(rec, source, remove_artifacts=False, lpf=True, bsub=True, lowpass=1000):
    """Perform a standardized strength analysis on a record selected by response_query or baseline_query.

    1. Determine timing of presynaptic stimulus pulse edges and spike
    2. Measure peak deflection on raw trace
    3. Apply deconvolution / artifact removal / lpf
    4. Measure peak deflection on deconvolved trace
    """
    data = Trace(rec.data, sample_rate=db.default_sample_rate)
    if source == 'pulse_response':
        # Find stimulus pulse edges for artifact removal
        start = rec.pulse_start - rec.rec_start
        pulse_times = [start, start + rec.pulse_dur]
        if rec.spike_time is None:
            # these pulses failed QC, but we analyze them anyway to make all data visible
            spike_time = 11e-3
        else:
            spike_time = rec.spike_time - rec.rec_start
    elif source == 'baseline':
        # Fake stimulus information to ensure that background data receives
        # the same filtering / windowing treatment
        pulse_times = [10e-3, 12e-3]
        spike_time = 11e-3
    else:
        raise ValueError("Invalid source %s" % source)

    results = {}

    results['raw_trace'] = data
    results['pulse_times'] = pulse_times
    results['spike_time'] = spike_time

    # Measure crosstalk from pulse onset
    p1 = data.time_slice(pulse_times[0]-200e-6, pulse_times[0]).median()
    p2 = data.time_slice(pulse_times[0], pulse_times[0]+200e-6).median()
    results['crosstalk'] = p2 - p1

    # crosstalk artifacts in VC are removed before deconvolution
    if rec.clamp_mode == 'vc' and remove_artifacts is True:
        data = remove_crosstalk_artifacts(data, pulse_times)
        remove_artifacts = False

    # Measure deflection on raw data
    results['pos_amp'], _ = measure_peak(data, '+', spike_time, pulse_times)
    results['neg_amp'], _ = measure_peak(data, '-', spike_time, pulse_times)

    # Deconvolution / artifact removal / filtering
    tau = 15e-3 if rec.clamp_mode == 'ic' else 5e-3
    dec_data = deconv_filter(data, pulse_times, tau=tau, lpf=lpf, remove_artifacts=remove_artifacts, bsub=bsub, lowpass=lowpass)
    results['dec_trace'] = dec_data

    # Measure deflection on deconvolved data
    results['pos_dec_amp'], results['pos_dec_latency'] = measure_peak(dec_data, '+', spike_time, pulse_times)
    results['neg_dec_amp'], results['neg_dec_latency'] = measure_peak(dec_data, '-', spike_time, pulse_times)
    
    return results
