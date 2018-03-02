# coding: utf8
"""
Big question: what's the best way to measure synaptic strength / connectivity?

1. Whatever method we use to measure connectivity, we also need to characterize the detection limit (per synapse)
2. Any method we choose should be run on both pulse response and background data, and the distributions of
   these results must be compared to make the connectivity call

"""
from __future__ import print_function, division

from collections import OrderedDict
import argparse, time, sys, os, pickle, io, multiprocessing
import numpy as np
import scipy.stats
import pandas

from sqlalchemy.orm import aliased

import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore

from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import Trace, TraceList
from neuroanalysis import filter
from neuroanalysis.event_detection import exp_deconvolve

from multipatch_analysis.database import database as db
from multipatch_analysis import config, synphys_cache
from multipatch_analysis.ui.multipatch_nwb_viewer import MultipatchNwbViewer
from multipatch_analysis.constants import EXCITATORY_CRE_TYPES, INHIBITORY_CRE_TYPES
import multipatch_analysis.qc as qc 



class TableGroup(object):
    def __init__(self):
        self.mappings = {}
        self.create_mappings()

    def __getitem__(self, item):
        return self.mappings[item]

    def create_mappings(self):
        for k,schema in self.schemas.items():
            self.mappings[k] = db.generate_mapping(k, schema)

    def drop_tables(self):
        for k in self.schemas:
            if k in db.engine.table_names():
                self[k].__table__.drop(bind=db.engine)

    def create_tables(self):
        for k in self.schemas:
            if k not in db.engine.table_names():
                self[k].__table__.create(bind=db.engine)


class PulseResponseStrengthTableGroup(TableGroup):
    """Measures pulse amplitudes for each pulse response and background chunk.
    """
    schemas = {
        'pulse_response_strength': [
            ('pulse_response_id', 'pulse_response.id', '', {'index': True}),
            ('pos_amp', 'float'),
            ('neg_amp', 'float'),
            ('pos_dec_amp', 'float'),
            ('neg_dec_amp', 'float'),
            ('pos_dec_latency', 'float'),
            ('neg_dec_latency', 'float'),
        ],
        'baseline_response_strength' : [
            ('baseline_id', 'baseline.id', '', {'index': True}),
            ('pos_amp', 'float'),
            ('neg_amp', 'float'),
            ('pos_dec_amp', 'float'),
            ('neg_dec_amp', 'float'),
            ('pos_dec_latency', 'float'),
            ('neg_dec_latency', 'float'),
        ]
        #'deconv_pulse_response': [
            #"Exponentially deconvolved pulse responses",
        #],
    }

    def create_mappings(self):
        TableGroup.create_mappings(self)
        
        PulseResponseStrength = self['pulse_response_strength']
        BaselineResponseStrength = self['baseline_response_strength']
        
        db.PulseResponse.pulse_response_strength = db.relationship(PulseResponseStrength, back_populates="pulse_response", cascade="delete", single_parent=True)
        PulseResponseStrength.pulse_response = db.relationship(db.PulseResponse, back_populates="pulse_response_strength", single_parent=True)

        db.Baseline.baseline_response_strength = db.relationship(BaselineResponseStrength, back_populates="baseline", cascade="delete", single_parent=True)
        BaselineResponseStrength.baseline = db.relationship(db.Baseline, back_populates="baseline_response_strength", single_parent=True)


class ConnectionStrengthTableGroup(TableGroup):
    schemas = {
        'connection_strength': [
            ('pair_id', 'pair.id', '', {'index': True}),
            ('synapse_type', 'str', '"ex" or "in"'),

            # current clamp metrics
            ('ic_n_samples', 'int'),
            # amplitude,
            ('ic_amp_med', 'float'),
            ('ic_amp_stdev', 'float'),
            ('ic_base_amp_med', 'float'),
            ('ic_base_amp_stdev', 'float'),
            ('ic_amp_ttest', 'float'),
            ('ic_amp_ks2samp', 'float'),
            # deconvolved amplitide
            ('ic_deconv_amp_med', 'float'),
            ('ic_deconv_amp_stdev', 'float'),
            ('ic_base_deconv_amp_med', 'float'),
            ('ic_base_deconv_amp_stdev', 'float'),
            ('ic_deconv_amp_ttest', 'float'),
            ('ic_deconv_amp_ks2samp', 'float'),
            # latency
            ('ic_latency_med', 'float'),
            ('ic_latency_stdev', 'float'),
            ('ic_base_latency_med', 'float'),
            ('ic_base_latency_stdev', 'float'),
            ('ic_latency_ttest', 'float'),
            ('ic_latency_ks2samp', 'float'),
            
            # voltage clamp metrics
            ('vc_n_samples', 'int'),
            # amplitude,
            ('vc_amp_med', 'float'),
            ('vc_amp_stdev', 'float'),
            ('vc_base_amp_med', 'float'),
            ('vc_base_amp_stdev', 'float'),
            ('vc_amp_ttest', 'float'),
            ('vc_amp_ks2samp', 'float'),
            # deconvolved amplitide
            ('vc_deconv_amp_med', 'float'),
            ('vc_deconv_amp_stdev', 'float'),
            ('vc_base_deconv_amp_med', 'float'),
            ('vc_base_deconv_amp_stdev', 'float'),
            ('vc_deconv_amp_ttest', 'float'),
            ('vc_deconv_amp_ks2samp', 'float'),
            # latency
            ('vc_latency_med', 'float'),
            ('vc_latency_stdev', 'float'),
            ('vc_base_latency_med', 'float'),
            ('vc_base_latency_stdev', 'float'),
            ('vc_latency_ttest', 'float'),
            ('vc_latency_ks2samp', 'float'),
            
        ],
    }

    def create_mappings(self):
        TableGroup.create_mappings(self)
        
        ConnectionStrength = self['connection_strength']
        
        db.Pair.connection_strength = db.relationship(ConnectionStrength, back_populates="pair", cascade="delete", single_parent=True)
        ConnectionStrength.pair = db.relationship(db.Pair, back_populates="connection_strength", single_parent=True)


pulse_response_strength_tables = PulseResponseStrengthTableGroup()
connection_strength_tables = ConnectionStrengthTableGroup()

def init_tables():
    global PulseResponseStrength, BaselineResponseStrength, ConnectionStrength
    pulse_response_strength_tables.create_tables()
    connection_strength_tables.create_tables()

    PulseResponseStrength = pulse_response_strength_tables['pulse_response_strength']
    BaselineResponseStrength = pulse_response_strength_tables['baseline_response_strength']
    ConnectionStrength = connection_strength_tables['connection_strength']


def measure_peak(trace, sign, spike_time, pulse_times, spike_delay=1e-3, response_window=4e-3):
    # Start measuring response after the pulse has finished, and no earlier than 1 ms after spike onset
    # response_start = max(spike_time + spike_delay, pulse_times[1])

    # Start measuring after spike and hope that the pulse offset doesn't get in the way
    # (if we wait for the pulse to end, then we miss too many fast rise / short latency events)
    response_start = max(spike_time + spike_delay, pulse_times[1])
    response_stop = response_start + response_window

    # measure baseline from beginning of data until 50µs before pulse onset
    baseline_start = 0
    baseline_stop = pulse_times[0] - 50e-6

    baseline = trace.time_slice(baseline_start, baseline_stop).data.mean()
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
        dt = dec.dt
        r = [-50e-6, 250e-6]
        edges = [(int((t+r[0])/dt), int((t+r[1])/dt)) for t in pulse_times]
        cleaned = filter.remove_artifacts(dec, edges, window=400e-6)
    else:
        cleaned = dec

    if bsub:
        baseline = np.median(cleaned.time_slice(cleaned.t0, cleaned.t0+10e-3).data)
        b_subbed = cleaned-baseline
    else:
        b_subbed = cleaned

    if lpf:
        return filter.bessel_filter(b_subbed, lowpass)
    else:
        return b_subbed


@db.default_session
def rebuild_strength(parallel=True, workers=6, session=None):
    for source in ['pulse_response', 'baseline']:
        print("Rebuilding %s strength table.." % source)
        
        # Divide workload into equal parts for each worker
        max_pulse_id = session.execute('select max(id) from %s' % source).fetchone()[0]
        chunk = 1 + (max_pulse_id // workers)
        parts = [(source, chunk*i, chunk*(i+1)) for i in range(workers)]

        if parallel:
            pool = multiprocessing.Pool(processes=workers)
            pool.map(compute_strength, parts)
        else:
            for part in parts:
                compute_strength(part)


def compute_strength(inds, session=None):
    # Thin wrapper just to allow calling from pool.map
    return _compute_strength(inds, session=session)


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

    # Ignore anything that failed QC
    q = q.filter(((db.PulseResponse.ex_qc_pass==True) | (db.PulseResponse.in_qc_pass==True)))

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

    # Ignore anything that failed QC
    q = q.filter(((db.Baseline.ex_qc_pass==True) | (db.Baseline.in_qc_pass==True)))

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


@db.default_session
def _compute_strength(inds, session=None):
    """Comput per-pulse-response strength metrics
    """
    source, start_id, stop_id = inds
    if source == 'baseline':
        q = baseline_query(session)
    elif source == 'pulse_response':
        q = response_query(session)
    else:
        raise ValueError("Invalid source %s" % source)

    prof = pg.debug.Profiler(delayed=False)
    
    next_id = start_id
    while True:
        # Request just a chunk of all pulse responses based on ID range
        if source == 'pulse_response':
            q1 = q.filter(db.PulseResponse.id>=next_id).filter(db.PulseResponse.id<stop_id)
            q1 = q1.order_by(db.PulseResponse.id)
        else:
            q1 = q.filter(db.Baseline.id>=next_id).filter(db.Baseline.id<stop_id)
            q1 = q1.order_by(db.Baseline.id)
        q1 = q1.limit(1000)  # process in 1000-record chunks

        prof('exec')
        recs = q1.all()
        prof('fetch')
        if len(recs) == 0:
            break
        new_recs = []

        for rec in recs:
            result = analyze_response_strength(rec, source)
            new_rec = {'%s_id'%source: rec.response_id}
            # copy a subset of results over to new record
            for k in ['pos_amp', 'neg_amp', 'pos_dec_amp', 'neg_dec_amp', 'pos_dec_latency','neg_dec_latency']:
                new_rec[k] = result[k]
            new_recs.append(new_rec)
        
        next_id = rec.response_id + 1
        
        sys.stdout.write("%d / %d\r" % (next_id-start_id, stop_id-start_id))
        sys.stdout.flush()

        prof('process')
        if source == 'pulse_response':
            session.bulk_insert_mappings(PulseResponseStrength, new_recs)
        else:
            session.bulk_insert_mappings(BaselineResponseStrength, new_recs)
        prof('insert')
        new_recs = []

        session.commit()
        prof('commit')


@db.default_session
def rebuild_connectivity(session):
    print("Rebuilding connectivity table..")
    
    expts_in_db = list_experiments(session=session)
    for i,expt in enumerate(expts_in_db):
        for pair in expt.pairs:
            devs = pair.pre_cell.electrode.device_id, pair.post_cell.electrode.device_id
            # Whether the user marked this as connected
            cell_ids = (devs[0] + 1, devs[1] + 1)
            
            fields = {}
            
            for clamp_mode in ['ic', 'vc']:
                amps = get_amps(session, pair, clamp_mode=clamp_mode)
                base_amps = get_baseline_amps(session, pair, limit=len(amps), clamp_mode=clamp_mode)
                if len(amps) == 0 or len(base_amps) == 0:
                    continue
                
                # decide whether to treat this connection as excitatory or inhibitory
                # (probably we can do much better here)
                pos_amp = amps['pos_dec_amp'].mean() - base_amps['pos_dec_amp'].mean()
                neg_amp = amps['neg_dec_amp'].mean() - base_amps['neg_dec_amp'].mean()
                if pos_amp > -neg_amp:
                    conn.synapse_type = 'ex'
                    pfx = 'pos_'
                    qc_field = 'ex_qc_pass' if clamp_mode == 'ic' else 'in_qc_pass'
                else:
                    conn.synapse_type = 'in'
                    pfx = 'neg_'
                    qc_field = 'ex_qc_pass' if clamp_mode == 'vc' else 'in_qc_pass'

                # filter out responses that failed qc
                amps_qc = amps[amps[qc_field]]
                base_qc = base_amps[base_amps[qc_field]]
                
                # select out positive or negative amplitude columns
                amp, base_amp = amps_qc[pfx+'amp'], base_qc[pfx+'amp']
                dec_amp, dec_base_amp = amps_qc[pfx+'dec_amp'], base_qc[pfx+'dec_amp']
                latency, base_latency = amps_qc[pfx+'dec_latency'], base_qc[pfx+'dec_latency']
        
                # compute mean/stdev of samples
                n_samp = len(amp)
                fields[clamp_mode + '_n_samples'] = n_samp
                if n_samp == 0:
                    continue
                
                # measure median, stdev, and statistical differences between
                # fg and bg for each measurement
                for (val, fg, bg) in [('amp', amp, base_amp), 
                                      ('deconv_amp', dec_amp, dec_base_amp),
                                      ('latency', latency, base_latency)]:
                    fields[clamp_mode + '_' + val + '_med'] = np.median(fg)
                    fields[clamp_mode + '_' + val + '_stdev'] = np.std(fg)
                    fields[clamp_mode + '_base_' + val + '_med'] = np.median(bg)
                    fields[clamp_mode + '_base_' + val + '_stdev'] = np.std(bg)
                    fields[clamp_mode + '_' + val + '_ttest'] = scipy.stats.ttest_ind(fg, bg, equal_var=False).pvalue
                    fields[clamp_mode + '_' + val + '_ks2samp'] = scipy.stats.ks_2samp(fg, bg).pvalue

            conn = ConnectionStrength(pair_id=pair.id, **fields)
            session.add(conn)
        
        session.commit()
        sys.stdout.write("%d / %d       \r" % (i, len(expts_in_db)))
        sys.stdout.flush()


@db.default_session
def list_experiments(session):
    return session.query(db.Experiment).all()


# @db.default_session
# def get_experiment_devs(expt, session):
#     devs = session.query(db.Recording.electr_id).join(db.SyncRec).filter(db.SyncRec.experiment_id==expt.id)
    
#     # Only keep channels with >2 recordings
#     counts = {}
#     for r in devs:
#         counts.setdefault(r[0], 0)
#         counts[r[0]] += 1
#     devs = [k for k,v in counts.items() if v > 2]
#     devs.sort()
    
#     return devs


# def get_experiment_pairs(expt):
#     devs = get_experiment_devs(expt)
#     return [(pre, post) for pre in devs for post in devs if pre != post]


class ExperimentBrowser(pg.TreeWidget):
    def __init__(self):
        pg.TreeWidget.__init__(self)
        self.setColumnCount(7)
        self.setHeaderLabels(['date', 'rig', 'organism', 'region', 'genotype', 'acsf'])
        self.populate()
        self._last_expanded = None
        
    def populate(self):
        self.items_by_pair_id = {}
        
        self.session = db.Session()
        db_expts = list_experiments(session=self.session)
        db_expts.sort(key=lambda e: e.acq_timestamp)
        for expt in db_expts:
            date = expt.acq_timestamp
            date_str = date.strftime('%Y-%m-%d')
            slice = expt.slice
            expt_item = pg.TreeWidgetItem(map(str, [date_str, expt.rig_name, slice.species, expt.target_region, slice.genotype, expt.acsf]))
            expt_item.expt = expt
            self.addTopLevelItem(expt_item)

            for pair in expt.pairs:
                if pair.n_ex_test_spikes == 0 and pair.n_in_test_spikes == 0:
                    continue
                cells = '%d => %d' % (pair.pre_cell.ext_id, pair.post_cell.ext_id)
                conn = {True:"connected", False:"unconnected", None:"?"}[pair.synapse]
                types = 'L%s %s => L%s %s' % (pair.pre_cell.target_layer or "?", pair.pre_cell.cre_type, pair.post_cell.target_layer or "?", pair.post_cell.cre_type)
                pair_item = pg.TreeWidgetItem([cells, conn, types])
                expt_item.addChild(pair_item)
                pair_item.pair = pair
                pair_item.expt = expt
                self.items_by_pair_id[pair.id] = pair_item
                
    def select(self, pair_id):
        """Select a specific pair from the list
        """
        if self._last_expanded is not None:
            self._last_expanded.setExpanded(False)
        item = self.items_by_pair_id[pair_id]
        self.clearSelection()
        item.setSelected(True)
        parent = item.parent()
        if not parent.isExpanded():
            parent.setExpanded(True)
            self._last_expanded = parent
        else:
            self._last_expanded = None
        self.scrollToItem(item)


def get_amps(session, pair, clamp_mode='ic'):
    """Select records from pulse_response_strength table
    """
    q = session.query(
        PulseResponseStrength.id,
        PulseResponseStrength.pos_amp,
        PulseResponseStrength.neg_amp,
        PulseResponseStrength.pos_dec_amp,
        PulseResponseStrength.neg_dec_amp,
        PulseResponseStrength.pos_dec_latency,
        PulseResponseStrength.neg_dec_latency,
        db.PulseResponse.ex_qc_pass,
        db.PulseResponse.in_qc_pass,
        db.PatchClampRecording.clamp_mode,
    ).join(db.PulseResponse)
    
    q, pre_rec, post_rec = join_pulse_response_to_expt(q)
        
    filters = [
        (pre_rec.electrode==pair.pre_cell.electrode,),
        (post_rec.electrode==pair.post_cell.electrode,),
        (db.PatchClampRecording.clamp_mode==clamp_mode,),
        (db.PatchClampRecording.qc_pass==True,),
    ]
    for filter_args in filters:
        q = q.filter(*filter_args)
    
    df = pandas.read_sql_query(q.statement, q.session.bind)
    recs = df.to_records()
    return recs


def get_baseline_amps_NO_ORM_VERSION(session, expt, dev, clamp_mode='ic'):
    # Tested this against the ORM version below; no difference in performance.
    query = """
    select 
        brs.id, brs.pos_amp, brs.neg_amp, brs.pos_dec_amp, brs.neg_dec_amp
    from 
        baseline_response_strength brs
    join baseline on brs.baseline_id = baseline.id
    join recording on baseline.recording_id = recording.id
    join patch_clamp_recording on patch_clamp_recording.recording_id = recording.id
    join sync_rec on recording.sync_rec_id = sync_rec.id
    join experiment on sync_rec.experiment_id = experiment.id
    where 
        experiment.id={expt_id} and
        recording.device_key={dev_key} and
        patch_clamp_recording.clamp_mode='{clamp_mode}' and
        patch_clamp_recording.baseline_potential<=-50e-3 and
        patch_clamp_recording.baseline_current>-800e-12 and
        patch_clamp_recording.baseline_current<800e-12
    """.format(expt_id=expt.id, dev_key=dev, clamp_mode=clamp_mode)
    df = pandas.read_sql(query, session.bind)
    recs = df.to_records()
    return recs


def get_baseline_amps(session, pair, clamp_mode='ic', limit=None):
    """Select records from baseline_response_strength table
    """
    q = session.query(
        BaselineResponseStrength.id,
        BaselineResponseStrength.pos_amp,
        BaselineResponseStrength.neg_amp,
        BaselineResponseStrength.pos_dec_amp,
        BaselineResponseStrength.neg_dec_amp,
        BaselineResponseStrength.pos_dec_latency,
        BaselineResponseStrength.neg_dec_latency,
        db.Baseline.ex_qc_pass,
        db.Baseline.in_qc_pass,
        db.PatchClampRecording.clamp_mode,
    ).join(db.Baseline).join(db.Recording).join(db.PatchClampRecording).join(db.SyncRec).join(db.Experiment)
    
    filters = [
        (db.Recording.electrode==pair.post_cell.electrode,),
        (db.PatchClampRecording.clamp_mode==clamp_mode,),
        (db.PatchClampRecording.qc_pass==True,),
    ]
    for filter_args in filters:
        q = q.filter(*filter_args)
    
    if limit is not None:
        q = q.limit(limit)

    df = pandas.read_sql_query(q.statement, q.session.bind)
    recs = df.to_records()
    return recs


def join_pulse_response_to_expt(query):
    pre_rec = aliased(db.Recording)
    post_rec = aliased(db.Recording)
    joins = [
        (post_rec, db.PulseResponse.recording),
        (db.PatchClampRecording,),
        (db.SyncRec,),
        (db.Experiment,),
        (db.StimPulse, db.PulseResponse.stim_pulse),
        (pre_rec, db.StimPulse.recording),
    ]
    for join_args in joins:
        query = query.join(*join_args)

    return query, pre_rec, post_rec


class ResponseStrengthPlots(QtGui.QWidget):
    def __init__(self, session):
        QtGui.QWidget.__init__(self)
        self.session = session

        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)

        self.analyses = [('neg', 'ic'), ('pos', 'ic'), ('neg', 'vc'), ('pos', 'vc')]
        self.analyzers = []
        for col, analysis in enumerate(self.analyses):
            analyzer = ResponseStrengthAnalyzer(analysis, session)
            self.analyzers.append(analyzer)
            self.layout.addWidget(analyzer.widget, 0, col)
                
    def load_conn(self, pair):
        with pg.BusyCursor():

            for analyzer in self.analyzers:
                analyzer.load_conn(pair)


class ResponseStrengthAnalyzer(object):
    def __init__(self, analysis, session):
        self.analysis = analysis  # ('pos'|'neg', 'ic'|'vc')
        self.session = session

        self.widget = QtGui.QWidget()
        self.layout = QtGui.QGridLayout()
        self.widget.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)

        self.gl = pg.GraphicsLayoutWidget()
        self.layout.addWidget(self.gl, 0, 0)

        # histogram plots
        self.hist_plot = pg.PlotItem(title=' '.join(analysis)+":")
        self.gl.addItem(self.hist_plot, row=0, col=0)
        self.hist_plot.hideAxis('bottom')

        # event scatter plots
        self.scatter_plot = pg.PlotItem()
        self.gl.addItem(self.scatter_plot, row=1, col=0)
        self.scatter_plot.getAxis('left').setTicks([[(1.4, 'fg'), (0.4, 'bg')], []])
        self.scatter_plot.setFixedHeight(200)
        self.hist_plot.setXLink(self.scatter_plot)

        self.bg_scatter = pg.ScatterPlotItem(pen=None, symbol='o', symbolPen=None)
        self.fg_scatter = pg.ScatterPlotItem(pen=None, symbol='o', symbolPen=None)
        self.scatter_plot.addItem(self.bg_scatter)
        self.scatter_plot.addItem(self.fg_scatter)
        self.fg_scatter.sigClicked.connect(self.fg_scatter_clicked)
        self.bg_scatter.sigClicked.connect(self.bg_scatter_clicked)

        # event selector
        self.event_region = pg.LinearRegionItem()
        self.scatter_plot.addItem(self.event_region)
        self.event_region.setZValue(-10)
        self.event_region.sigRegionChangeFinished.connect(self.region_changed)

        # trace plots
        self.fg_trace_plot = pg.PlotItem()
        self.gl.addItem(self.fg_trace_plot, row=2, col=0)
        self.fg_trace_plot.hideAxis('bottom')
        self.fg_trace_plot.addLine(x=10e-3, pen=(0, 0, 100), movable=False)
        self.bg_trace_plot = pg.PlotItem(labels={'bottom': ('time', 's')})
        self.gl.addItem(self.bg_trace_plot, row=3, col=0)
        self.bg_trace_plot.setXLink(self.fg_trace_plot)
        self.bg_trace_plot.setYLink(self.fg_trace_plot)
        self.bg_trace_plot.addLine(x=10e-3, pen=(0, 0, 100), movable=False)
        self.fg_trace_plot.setXRange(0, 20e-3)

        self.ctrl = QtGui.QWidget()
        self.layout.addWidget(self.ctrl, 1, 0)
        self.ctrl_layout = QtGui.QGridLayout()
        self.ctrl.setLayout(self.ctrl_layout)
        self.ctrl_layout.setContentsMargins(0, 0, 0, 0)

        self.deconv_check = QtGui.QCheckBox('deconvolve')
        self.deconv_check.setChecked(True)
        self.ctrl_layout.addWidget(self.deconv_check, 0, 0)
        self.deconv_check.toggled.connect(self.replot_all)

        self.bsub_check = QtGui.QCheckBox('bsub')
        self.bsub_check.setChecked(True)
        self.ctrl_layout.addWidget(self.bsub_check, 0, 1)
        self.bsub_check.toggled.connect(self.replot_all)

        self.lpf_check = QtGui.QCheckBox('lpf')
        self.lpf_check.setChecked(True)
        self.ctrl_layout.addWidget(self.lpf_check, 0, 2)
        self.lpf_check.toggled.connect(self.replot_all)

        self.ar_check = QtGui.QCheckBox('crosstalk')
        self.ar_check.setChecked(True)
        self.ctrl_layout.addWidget(self.ar_check, 0, 3)
        self.ar_check.toggled.connect(self.replot_all)

        self.align_check = QtGui.QCheckBox('align')
        self.align_check.setChecked(True)
        self.ctrl_layout.addWidget(self.align_check, 0, 4)
        self.align_check.toggled.connect(self.replot_all)

        self.selected_fg_traces = []
        self.selected_bg_traces = []
        self.clicked_fg_traces = []
        self.clicked_bg_traces = []
        self._selected_fg_ids = []
        self._selected_bg_ids = []
        self._clicked_fg_ids = []
        self._clicked_bg_ids = []
    
    def fg_scatter_clicked(self, sp, points):
        """Point(s) were clicked; plot their source traces in a different color.
        """
        ids = [p.data()['id'] for p in points]
        self._clicked_fg_ids = ids
        self.plot_prd_ids(ids, 'fg', trace_list=self.clicked_fg_traces, pen='y')

        global selected_response
        selected_response = self.session.query(PulseResponseStrength).filter(PulseResponseStrength.id==ids[0]).first()
        prs_qc()

    def bg_scatter_clicked(self, sp, points):
        """Point(s) were clicked; plot their source traces in a different color.
        """
        ids = [p.data()['id'] for p in points]
        self._clicked_bg_ids = ids
        self.plot_prd_ids(ids, 'bg', trace_list=self.clicked_bg_traces, pen='y')

        global selected_response
        selected_response = self.session.query(BaselineResponseStrength).filter(BaselineResponseStrength.id==ids[0]).first()
        prs_qc()

    def load_data(self):
        amp_recs = get_amps(self.session, self.pair, clamp_mode=self.analysis[1])
        base_recs = get_baseline_amps(self.session, self.pair, limit=len(amp_recs), clamp_mode=self.analysis[1])
        return amp_recs, base_recs

    def load_conn(self, pair):
        self.pair = pair
        amp_recs, base_recs = self.load_data()
        
        # select fg/bg data
        fg_data = amp_recs
        bg_data = base_recs[:len(fg_data)]
        if self.analysis[0] == 'pos':
            data_field = 'pos_dec_amp'
            qc_field = 'ex_qc_pass' if self.analysis[1] == 'ic' else 'in_qc_pass'
        elif self.analysis[0] == 'neg':
            data_field = 'neg_dec_amp'
            qc_field = 'in_qc_pass' if self.analysis[1] == 'ic' else 'ex_qc_pass'
        fg_x = fg_data[data_field]
        bg_x = bg_data[data_field]
        fg_qc = fg_data[qc_field]
        bg_qc = bg_data[qc_field]

        # record data for later use
        self.fg_x = fg_x
        self.bg_x = bg_x
        self.fg_qc = fg_qc
        self.bg_qc = bg_qc
        self.qc_field = qc_field
        self.fg_data = fg_data
        self.bg_data = bg_data

        # select colors
        pass_brush = pg.mkBrush((255, 255, 255, 80))
        fail_brush = pg.mkBrush((150, 0, 0, 80))
        fg_color = [(pass_brush if qc_pass else fail_brush) for qc_pass in fg_qc]
        bg_color = [(pass_brush if qc_pass else fail_brush) for qc_pass in bg_qc]
        
        # clear old plots
        self.fg_scatter.setData([])
        self.bg_scatter.setData([])
        self.hist_plot.clear()
        self.hist_plot.setTitle("")
        self.clear_trace_plots()
        # clear click-selected plots
        self._clicked_fg_ids = []
        self._clicked_bg_ids = []

        if len(fg_x) == 0:
            return

        # scatter plots of fg/bg data
        fg_y = 1 + np.random.random(size=len(fg_x)) * 0.8
        bg_y = np.random.random(size=len(bg_x)) * 0.8
        self.fg_scatter.setData(fg_x, fg_y, data=fg_data, brush=fg_color)
        self.bg_scatter.setData(bg_x, bg_y, data=bg_data, brush=bg_color)

        # show only qc-passed data in histogram
        fg_x_qc = fg_x[fg_qc] 
        bg_x_qc = bg_x[bg_qc]
        if len(fg_x_qc) == 0 or len(bg_x_qc) == 0:
            return
        
        # plot histograms        
        n_bins = max(5, int(len(fg_x_qc)**0.5))
        all_x = np.concatenate([fg_x_qc, bg_x_qc])
        bins = np.linspace(np.percentile(all_x, 2), np.percentile(all_x, 98), n_bins+1)
        fg_hist = np.histogram(fg_x_qc, bins=bins)
        bg_hist = np.histogram(bg_x_qc, bins=bins)
        self.hist_plot.plot(bg_hist[1], bg_hist[0], stepMode=True, fillLevel=0, brush=(255, 0, 0, 100), pen=0.5)
        self.hist_plot.plot(fg_hist[1], fg_hist[0], stepMode=True, fillLevel=0, brush=(0, 255, 255, 100), pen=1.0)

        # measure distance between distributions
        ks = scipy.stats.ks_2samp(fg_x_qc, bg_x_qc)
        self.hist_plot.setTitle("%s: KS=%0.02g" % (' '.join(self.analysis), ks.pvalue))

        # set event region, but don't update (it's too expensive)
        try:
            self.event_region.sigRegionChangeFinished.disconnect(self.region_changed)
            self.event_region.setRegion([fg_x_qc.min(), fg_x_qc.max()])
        finally:
            self.event_region.sigRegionChangeFinished.connect(self.region_changed)

    def clear_trace_plots(self):
        self.fg_trace_plot.clear()
        self.bg_trace_plot.clear()
        self.selected_fg_traces = []
        self.selected_bg_traces = []

    def region_changed(self):
        """Event selection region changed; replot selected traces
        """
        self.clear_trace_plots()
        rgn = self.event_region.getRegion()

        # Find selected events
        mask = (self.fg_x > rgn[0]) & (self.fg_x < rgn[1])
        fg_ids = self.fg_data['id'][mask]
        self._selected_fg_ids = fg_ids

        mask = (self.bg_x > rgn[0]) & (self.bg_x < rgn[1])
        bg_ids = self.bg_data['id'][mask]
        self._selected_bg_ids = bg_ids

        # generate trace queries
        self.plot_prd_ids(fg_ids, 'fg', trace_list=self.selected_fg_traces, avg=True)
        self.plot_prd_ids(bg_ids, 'bg', trace_list=self.selected_bg_traces, avg=True)

    def replot_all(self):
        self.plot_prd_ids(self._selected_fg_ids, 'fg', pen=None, trace_list=self.selected_fg_traces, avg=True)
        self.plot_prd_ids(self._selected_bg_ids, 'bg', pen=None, trace_list=self.selected_bg_traces, avg=True)
        self.plot_prd_ids(self._clicked_fg_ids, 'fg', pen='y', trace_list=self.clicked_fg_traces)
        self.plot_prd_ids(self._clicked_bg_ids, 'bg', pen='y', trace_list=self.clicked_bg_traces)

    def plot_prd_ids(self, ids, source, pen=None, trace_list=None, avg=False):
        """Plot raw or decolvolved PulseResponse data, given IDs of records in
        a PulseResponseStrength table.
        """
        with pg.BusyCursor():
            if source == 'fg':
                q = response_query(self.session)
                q = q.join(PulseResponseStrength)
                q = q.filter(PulseResponseStrength.id.in_(ids))
                q = q.add_column(db.PulseResponse.start_time)
                traces = self.selected_fg_traces
                plot = self.fg_trace_plot
            else:
                q = baseline_query(self.session)
                q = q.join(BaselineResponseStrength)
                q = q.filter(BaselineResponseStrength.id.in_(ids))
                q = q.add_column(db.Baseline.start_time)
                traces = self.selected_bg_traces
                plot = self.bg_trace_plot
            
            q = q.join(db.SyncRec).add_column(db.SyncRec.ext_id.label('sync_rec_ext_id'))
            recs = q.all()
            if len(recs) == 0:
                return
            
            for i in trace_list[:]:
                plot.removeItem(i)
                trace_list.remove(i)
                
            if pen is None:
                alpha = np.clip(1000 / len(recs), 30, 255)
                pen = (255, 255, 255, alpha)
                
            traces = []
            spike_times = []
            spike_values = []
            for rec in recs:
                # Filter by QC unless we selected just a single record
                if len(recs) > 0 and getattr(rec, self.qc_field) is False:
                    continue

                s = {'fg': 'pulse_response', 'bg': 'baseline'}[source]
                result = analyze_response_strength(rec, source=s, lpf=self.lpf_check.isChecked(), 
                                                   remove_artifacts=self.ar_check.isChecked(), bsub=self.bsub_check.isChecked())

                if self.deconv_check.isChecked():
                    trace = result['dec_trace']
                else:
                    trace = result['raw_trace']
                    if self.bsub_check.isChecked():
                        trace = trace - np.median(trace.time_slice(0, 9e-3).data)
                    if self.lpf_check.isChecked():
                        trace = filter.bessel_filter(trace, 500)
                
                spike_values.append(trace.value_at([result['spike_time']])[0])
                if self.align_check.isChecked():
                    trace.t0 = -result['spike_time']
                    spike_times.append(0)
                else:
                    spike_times.append(result['spike_time'])

                traces.append(trace)
                trace_list.append(plot.plot(trace.time_values, trace.data, pen=pen))

            if avg:
                mean = TraceList(traces).mean()
                trace_list.append(plot.plot(mean.time_values, mean.data, pen='g'))
                trace_list[-1].setZValue(10)

            spike_scatter = pg.ScatterPlotItem(spike_times, spike_values, size=4, pen=None, brush=(200, 200, 0))
            spike_scatter.setZValue(-100)
            plot.addItem(spike_scatter)
            trace_list.append(spike_scatter)


def query_all_pairs():
    query = ("""
    select """
        # ((DATE_PART('day', experiment.acq_timestamp - '1970-01-01'::timestamp) * 24 + 
        # DATE_PART('hour', experiment.acq_timestamp - '1970-01-01'::timestamp)) * 60 +
        # DATE_PART('minute', experiment.acq_timestamp - '1970-01-01'::timestamp)) * 60 +
        # DATE_PART('second', experiment.acq_timestamp - '1970-01-01'::timestamp) as acq_timestamp,
        """
        connection_strength.*,
        experiment.id as experiment_id,
        experiment.acq_timestamp as acq_timestamp,
        experiment.rig_name,
        experiment.acsf,
        slice.species,
        slice.genotype,
        slice.age,
        slice.slice_time,
        pre_cell.ext_id as pre_cell_id,
        pre_cell.cre_type as pre_cre_type,
        pre_cell.target_layer as pre_target_layer,
        post_cell.ext_id as post_cell_id,
        post_cell.cre_type as post_cre_type,
        post_cell.target_layer as post_target_layer,
        pair.synapse,
        pair.crosstalk_artifact,
        abs(post_cell.ext_id - pre_cell.ext_id) as electrode_distance
    from connection_strength
    join pair on connection_strength.pair_id=pair.id
    join cell pre_cell on pair.pre_cell_id=pre_cell.id
    join cell post_cell on pair.post_cell_id=post_cell.id
    join experiment on pair.expt_id=experiment.id
    join slice on experiment.slice_id=slice.id
    order by acq_timestamp
    """)
    session = db.Session()
    df = pandas.read_sql(query, session.bind)

    # test out a new metric:
    # df['connection_signal'] = pandas.Series(df['deconv_amp_med'] / df['latency_stdev'], index=df.index)
    # df['connection_background'] = pandas.Series(df['deconv_base_amp_med'] / df['base_latency_stdev'], index=df.index)

    ts = [datetime_to_timestamp(t) for t in df['acq_timestamp']]
    df['acq_timestamp'] = ts
    recs = df.to_records()
    return recs


def datetime_to_timestamp(d):
    return time.mktime(d.timetuple()) + d.microsecond * 1e-6


def prs_qc():
    """Convenience function for debugging QC: returns (recording, window) arguments
    used for pulse response QC
    """
    global selected_response
    sr = selected_response
    if isinstance(sr, PulseResponseStrength):
        resp = sr.pulse_response
    else:
        resp = sr.baseline
    rec = resp.recording
    start = resp.start_time
    expt = rec.sync_rec.experiment
    nwb = expt.data
    rrec = nwb.contents[rec.sync_rec.ext_id][rec.electrode.device_id]
    dt = rrec['primary'].dt
    dur = len(resp.data) / db.default_sample_rate
    w = [int(start/dt), int((start+dur)/dt)]

    stdev = rrec[w[0]:w[1]]['primary'].std()
    if rrec.clamp_mode == 'ic':
        median = rrec[w[0]:w[1]]['primary'].median()
    else:
        median = rrec[w[0]:w[1]]['command'].median()
    print("Selected:  sweep %d @ %d ms   stdev=%0.3g, med=%0.3g" % (rec.sync_rec.ext_id, np.round(start*1000), stdev, median))

    return rrec, w


init_tables()



if __name__ == '__main__':
    #tt = pg.debug.ThreadTrace()
    parser = argparse.ArgumentParser()
    parser.add_argument('--rebuild', action='store_true', default=False)
    parser.add_argument('--rebuild-connectivity', action='store_true', default=False, dest='rebuild_connectivity')
    parser.add_argument('--local', action='store_true', default=False)
    parser.add_argument('--workers', type=int, default=6)
    
    args, extra = parser.parse_known_args(sys.argv[1:])


    from multipatch_analysis.ui.multipatch_nwb_viewer import MultipatchNwbViewer    
    from multipatch_analysis.experiment_list import cached_experiments
    expts = cached_experiments()

    pg.dbg()

    if args.rebuild:
        connection_strength_tables.drop_tables()
        pulse_response_strength_tables.drop_tables()
        init_tables()
        rebuild_strength(parallel=(not args.local), workers=args.workers)
        rebuild_connectivity()
    elif args.rebuild_connectivity:
        print("drop tables..")
        connection_strength_tables.drop_tables()
        print("create tables..")
        init_tables()
        print("rebuild..")
        rebuild_connectivity()
    else:
        init_tables()
        
    
    # global session for querying from DB
    session = db.Session()
    
    win = pg.QtGui.QSplitter()
    win.setOrientation(pg.QtCore.Qt.Horizontal)
    win.resize(1000, 800)
    win.show()
    
    b = ExperimentBrowser()
    win.addWidget(b)
    
    rs_plots = ResponseStrengthPlots(session)
    win.addWidget(rs_plots)
    
    def selected(*args):
        """A pair was selected; update the event plots
        """
        global sel
        sel = b.selectedItems()
        if len(sel) == 0:
            return
        sel = sel[0]
        if hasattr(sel, 'pair'):
            pair = sel.pair
            expt = pair.experiment
            rs_plots.load_conn(pair)
        print(sel.expt.original_path)
        ts = sel.expt.acq_timestamp
        sec = datetime_to_timestamp(ts)
        print(sec)

    b.itemSelectionChanged.connect(selected)

    nwb_viewer = MultipatchNwbViewer()

    def dbl_clicked(index):
        with pg.BusyCursor():
            item = b.itemFromIndex(index)[0]
            nwb = os.path.join(item.expt.original_path, item.expt.ephys_file)
            print(nwb)
            try:
                cache = synphys_cache.get_cache().get_cache(nwb)
                print("load cached:", cache)
                nwb_viewer.load_nwb(cache)
            except Exception:
                print("load remote:", nwb)
                nwb_viewer.load_nwb(nwb)
            nwb_viewer.show()
        
    b.doubleClicked.connect(dbl_clicked)

    spw = pg.ScatterPlotWidget()
    spw.style['symbolPen'] = None
    
    spw.show()

    recs = query_all_pairs()

    fields = [
        ('synapse', {'mode': 'enum', 'values': [True, False, None]}),
        ('synapse_type', {'mode': 'enum', 'values': ['in', 'ex']}),
        ('pre_cre_type', {'mode': 'enum', 'values': list(set(recs['pre_cre_type']))}),
        ('post_cre_type', {'mode': 'enum', 'values': list(set(recs['post_cre_type']))}),
        ('pre_target_layer', {'mode': 'enum'}),
        ('post_target_layer', {'mode': 'enum'}),
        ('n_samples', {}),
        ('rig_name', {'mode': 'enum', 'values': list(set(recs['rig_name']))}),
        ('acsf', {'mode': 'enum', 'values': list(set(recs['acsf']))}),
        ('acq_timestamp', {}),
        ('crosstalk_artifact', {'units': 'V'}),
        ('electrode_distance', {}),
    ]
    fnames = [f[0] for f in fields]
    for f in recs.dtype.names:
        if f in fnames:
            continue
        if 'amp' in f:
            fields.append((f, {'units': 'V'}))
        elif 'latency' in f:
            fields.append((f, {'units': 's'}))
        else:
            fields.append((f, {}))
    spw.setFields(fields)
    spw.setData(recs)

    def conn_clicked(spw, points):
        d = points[0].data()
        b.select(pair_id=d['pair_id'])

    spw.sigScatterPlotClicked.connect(conn_clicked)



    # attempt a little machine learning. 
    # this could work if we can generate better features:
    #     KS measured on latency
    #     Gaussian fit to deconvolved trace (amp, time, width, nrmse)
    #         Seed by fitting average first
    #     KS on fit parameters
    #     Maybe something better than KS? Mixed gaussian model?
    #          Histogram values on the above features?  (might need a lot more data for this)
    from sklearn import svm, preprocessing

    features = recs[[
        'ic_n_samples', 
        #'amp_med', 'amp_stdev', 'base_amp_med', 'base_amp_stdev', 'amp_med_minus_base', 'amp_stdev_minus_base', 
        #'deconv_amp_med', 'deconv_amp_stdev', 'deconv_base_amp_med', 'deconv_base_amp_stdev', 'deconv_amp_med_minus_base', 'deconv_amp_stdev_minus_base', 
        #'amp_ttest',
        #'deconv_amp_ttest',
        #'amp_ks2samp', 
        'ic_deconv_amp_ks2samp',
        'ic_latency_med', 'ic_latency_stdev', 
        'ic_base_latency_med', 'ic_base_latency_stdev',
        #'abs_amp_med', 'abs_base_amp_med', 'abs_amp_med_minus_base', 
        'ic_deconv_amp_med', 'ic_base_deconv_amp_med', #'abs_deconv_amp_med_minus_base', 
        'electrode_distance'
    ]]

    features['ic_deconv_amp_ks2samp'] = np.log(-np.log(features['ic_deconv_amp_ks2samp']))
    # for k in ['deconv_amp_ks2samp']:
    #     features[k] = np.log(features[k])

    mask = features['ic_n_samples'] > 100

    x = np.array([tuple(r) for r in features[mask]])
    y = recs['synapse'][mask]

    order = np.arange(len(y))
    np.random.shuffle(order)
    x = x[order]
    y = y[order]

    mask = np.all(np.isfinite(x), axis=1) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    x = preprocessing.normalize(x, axis=0)

    # train on equal parts connected and non-connected
    train_mask = np.zeros(len(y), dtype='bool')
    syn = np.argwhere(y)
    n = len(syn) // 4 * 3
    train_mask[syn[:n]] = True
    other = np.arange(len(y))[~train_mask]
    np.random.shuffle(other)
    train_mask[other[:n]] = True

    train_x = x[train_mask]
    train_y = y[train_mask]
    test_x = x[~train_mask]
    test_y = y[~train_mask]
    print("Train: %d  test: %d" % (len(train_y), len(test_y)))

    clf = svm.LinearSVC()
    clf.fit(train_x, train_y)

    def test(clf, test_x, test_y):
        pred = clf.predict(test_x)
        def pr(name, pred, test_y):
            print("%s:  %d/%d  %0.2f%%" % (name, (pred & test_y).sum(), test_y.sum(), 100 * (pred & test_y).sum() / test_y.sum()))
            
        pr(" True positive", pred, test_y)
        pr("False positive", pred, ~test_y)
        pr(" True negative", ~pred, ~test_y)
        pr("False negative", ~pred, test_y)

    print("TRAIN:")
    test(clf, train_x, train_y)

    print("TEST:")
    test(clf, test_x, test_y)

    