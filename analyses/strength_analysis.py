"""
Big question: what's the best way to measure synaptic strength / connectivity?

"""
from __future__ import print_function, division

from collections import OrderedDict
from multipatch_analysis.constants import EXCITATORY_CRE_TYPES, INHIBITORY_CRE_TYPES

import argparse, time, sys, os, pickle, io, multiprocessing
import numpy as np
import scipy.stats

from sqlalchemy.orm import aliased

import pyqtgraph as pg

from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.filter import bessel_filter
from neuroanalysis.event_detection import exp_deconvolve

from multipatch_analysis.database import database as db
from multipatch_analysis import config
from multipatch_analysis.ui.multipatch_nwb_viewer import MultipatchNwbViewer


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


class ConnectionStrengthTableGroup(TableGroup):
    schemas = {
        'connection_strength': [
            ('pair_id', 'pair.id', '', {'index': True}),
            ('synapse_type', 'str', '"ex" or "in"'),
            ('n_samples', 'int'),
            ('amp_med', 'float'),
            ('amp_stdev', 'float'),
            ('base_amp_med', 'float'),
            ('base_amp_stdev', 'float'),
            ('amp_med_minus_base', 'float'),
            ('amp_stdev_minus_base', 'float'),
            ('deconv_amp_med', 'float'),
            ('deconv_amp_stdev', 'float'),
            ('deconv_base_amp_med', 'float'),
            ('deconv_base_amp_stdev', 'float'),
            ('deconv_amp_med_minus_base', 'float'),
            ('deconv_amp_stdev_minus_base', 'float'),
            ('amp_ttest', 'float'),
            ('deconv_amp_ttest', 'float'),
            ('amp_ks2samp', 'float'),
            ('deconv_amp_ks2samp', 'float'),
            ('latency_med', 'float'),
            ('latency_stdev', 'float'),
            ('base_latency_med', 'float'),
            ('base_latency_stdev', 'float'),
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


def measure_peak(trace, sign, baseline=(0e-3, 9e-3), response=(11e-3, 17e-3)):
    baseline = trace.time_slice(*baseline).data.mean()
    response = trace.time_slice(*response)
    if sign == '+':
        i = np.argmax(response.data)
    else:
        i = np.argmin(response.data)
    peak = response.data[i]
    latency = response.time_values[i]
    return peak - baseline, latency


def measure_sum(trace, sign, baseline=(0e-3, 9e-3), response=(12e-3, 17e-3)):
    baseline = trace.time_slice(*baseline).data.sum()
    peak = trace.time_slice(*response).data.sum()
    return peak - baseline

        
def deconv_filter(trace, tau=15e-3, lowpass=300.):
    dec = exp_deconvolve(trace, tau)
    baseline = np.median(dec.time_slice(dec.t0, dec.t0+10e-3).data)
    deconv = bessel_filter(dec-baseline, lowpass)
    return deconv


@db.default_session
def rebuild_strength(parallel=True, workers=6, session=None):
    for source in ['baseline', 'pulse_response']:
        print("Rebuilding %s strength table.." % source)
        
        max_pulse_id = session.execute('select max(id) from %s' % source).fetchone()[0]
        chunk = 1 + (max_pulse_id // workers)
        parts = [(source, chunk*i, chunk*(i+1)) for i in range(4)]

        if parallel:
            pool = multiprocessing.Pool(processes=workers)
            pool.map(compute_strength, parts)
        else:
            for part in parts:
                compute_strength(part)


def compute_strength(inds, session=None):
    # Thin wrapper just to allow calling from pool.map
    return _compute_strength(inds, session=session)


@db.default_session
def _compute_strength(inds, session=None):
    source, start_id, stop_id = inds
    q = """
        SELECT id, data FROM %s
        WHERE id >= %d and id < %d
        ORDER BY id
        LIMIT 1000
    """
    
    prof = pg.debug.Profiler(delayed=False)
    
    next_id = start_id
    while True:
        response = session.execute(q % (source, next_id, stop_id))  # bottleneck here ~30 ms
        prof('exec')
        recs = response.fetchall()
        prof('fetch')
        if len(recs) == 0:
            break
        new_recs = []
        for rec in recs:
            pr_id, data = rec
            data = np.load(io.BytesIO(data))
            #prof('load')
            
            new_rec = {'%s_id'%source: pr_id}
            
            data = Trace(data, sample_rate=20e3)
            new_rec['pos_amp'], _ = measure_peak(data, '+')
            new_rec['neg_amp'], _ = measure_peak(data, '-')
            #prof('comp 1')
            
            dec_data = deconv_filter(data)
            new_rec['pos_dec_amp'], new_rec['pos_dec_latency'] = measure_peak(dec_data, '+')
            new_rec['neg_dec_amp'], new_rec['neg_dec_latency'] = measure_peak(dec_data, '-')
            #prof('comp 3')
            
            new_recs.append(new_rec)
        
        next_id = pr_id + 1
        
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
            amps = get_amps(session, expt, pair)
            base_amps = get_baseline_amps(session, expt, pair, limit=len(amps))
            
            conn = ConnectionStrength(pair_id=pair.id)
            # Whether the user marked this as connected
            cell_ids = (devs[0] + 1, devs[1] + 1)
            
            # decide whether to treat this connection as excitatory or inhibitory
            # (probably we can do much better here)
            pos_amp = amps['pos_dec_amp'].mean() - base_amps['pos_dec_amp'].mean()
            neg_amp = amps['neg_dec_amp'].mean() - base_amps['neg_dec_amp'].mean()
            if pos_amp > -neg_amp:
                conn.synapse_type = 'ex'
                pfx = 'pos_'
            else:
                conn.synapse_type = 'in'
                pfx = 'neg_'
            # select out positive or negative amplitude columns
            amp, base_amp = amps[pfx+'amp'], base_amps[pfx+'amp']
            dec_amp, dec_base_amp = amps[pfx+'dec_amp'], base_amps[pfx+'dec_amp']
            latency, base_latency = amps[pfx+'dec_latency'], base_amps[pfx+'dec_latency']
    
            # compute mean/stdev of samples
            n_samp = len(amps)
            conn.n_samples = n_samp
            if n_samp == 0:
                continue
            conn.amp_med = np.median(amp)
            conn.amp_stdev = amp.std()
            conn.base_amp_med = np.median(base_amp)
            conn.base_amp_stdev = base_amp.std()
            conn.amp_med_minus_base = conn.amp_med - conn.base_amp_med
            conn.amp_stdev_minus_base = conn.amp_stdev - conn.base_amp_stdev
            conn.deconv_amp_med = np.median(dec_amp)
            conn.deconv_amp_stdev = dec_amp.std()
            conn.deconv_base_amp_med = np.median(dec_base_amp)
            conn.deconv_base_amp_stdev = dec_base_amp.std()
            conn.deconv_amp_med_minus_base = conn.deconv_amp_med - conn.deconv_base_amp_med
            conn.deconv_amp_stdev_minus_base = conn.deconv_amp_stdev - conn.deconv_base_amp_stdev

            # do some statistical tests
            conn.amp_ks2samp = scipy.stats.ks_2samp(amp, base_amp).pvalue
            conn.deconv_amp_ks2samp = scipy.stats.ks_2samp(dec_amp, dec_base_amp).pvalue
            conn.amp_ttest = scipy.stats.ttest_ind(amp, base_amp, equal_var=False).pvalue
            conn.deconv_amp_ttest = scipy.stats.ttest_ind(dec_amp, dec_base_amp, equal_var=False).pvalue

            # deconvolved peak latencies
            conn.latency_med = np.median(latency)
            conn.latency_stdev = np.std(latency)
            conn.base_latency_med = np.median(base_latency)
            conn.base_latency_stdev = np.std(base_latency)

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
        self.setColumnCount(6)
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
            expt_item = pg.TreeWidgetItem(map(str, [date_str, expt.rig_name, slice.species, expt.target_region, slice.genotype]))
            expt_item.expt = expt
            self.addTopLevelItem(expt_item)

            for pair in expt.pairs:
                pair_item = pg.TreeWidgetItem(['%d => %d' % (pair.pre_cell.ext_id, pair.post_cell.ext_id), str(pair.synapse)])
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
    ).join(db.PulseResponse)

    dtype = [
        ('id', 'int'),
        ('pos_amp', 'float'),
        ('neg_amp', 'float'),
        ('pos_dec_amp', 'float'),
        ('neg_dec_amp', 'float'),
        ('pos_dec_latency', 'float'),
        ('neg_dec_latency', 'float'),
    ]
    
    q, pre_rec, post_rec = join_pulse_response_to_expt(q)
        
    filters = [
        (pre_rec.electrode==pair.pre_cell.electrode,),
        (post_rec.electrode==pair.post_cell.electrode,),
        (db.PatchClampRecording.clamp_mode==clamp_mode,),
        (db.PatchClampRecording.baseline_potential<=-50e-3,),
        (db.PatchClampRecording.baseline_current>-800e-12,),
        (db.PatchClampRecording.baseline_current<400e-12,),
    ]
    for filter_args in filters:
        q = q.filter(*filter_args)
    
    # fetch all
    recs = q.all()
    
    # fill numpy array
    arr = np.empty(len(recs), dtype=dtype)
    for i,rec in enumerate(recs):
        arr[i] = rec

    return arr


def get_baseline_amps_NO_ORM_VERSION(session, expt, dev, clamp_mode='ic'):
    # Tested this against the ORM version below; no difference in performance.
    import pandas
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
    ).join(db.Baseline).join(db.Recording).join(db.PatchClampRecording).join(db.SyncRec).join(db.Experiment)

    dtype = [
        ('id', 'int'),
        ('pos_amp', 'float'),
        ('neg_amp', 'float'),
        ('pos_dec_amp', 'float'),
        ('neg_dec_amp', 'float'),
        ('pos_dec_latency', 'float'),
        ('neg_dec_latency', 'float'),
    ]
    
    filters = [
        (db.Recording.electrode==pair.post_cell.electrode,),
        (db.PatchClampRecording.clamp_mode==clamp_mode,),
        (db.PatchClampRecording.baseline_potential<=-50e-3,),
        (db.PatchClampRecording.baseline_current>-800e-12,),
        (db.PatchClampRecording.baseline_current<400e-12,),
    ]
    for filter_args in filters:
        q = q.filter(*filter_args)
    
    if limit is not None:
        q = q.limit(limit)

    # fetch all
    recs = q.all()
    
    # fill numpy array
    arr = np.empty(len(recs), dtype=dtype)
    for i,rec in enumerate(recs):
        arr[i] = rec

    return arr


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


class ResponseStrengthPlots(object):
    def __init__(self, session):
        self.session = session
        
        self.str_plot = pg.PlotItem()

        self.amp_sp = pg.ScatterPlotItem(pen=None, symbol='o', symbolPen=None)
        self.str_plot.addItem(self.amp_sp)
        self.amp_sp.sigClicked.connect(self.amp_clicked)

        self.base_sp = pg.ScatterPlotItem(pen=None, symbol='o', symbolPen=None)
        self.str_plot.addItem(self.base_sp)
        self.base_sp.sigClicked.connect(self.base_clicked)
        
        self.trace_plot = pg.PlotItem(labels={'bottom': ('time', 's'), 'left': ('Vm', 'V')})
        self.dec_plot = pg.PlotItem(labels={'bottom': ('time', 's'), 'left': ('Vm', 'V')})
    
    def amp_clicked(self, sp, points):
        self.clicked('amp', points)

    def base_clicked(self, sp, points):
        self.clicked('base', points)

    def clicked(self, source, points):
        self.clicked_points = points
        
        ids = [p.data()['id'] for p in points]
        if source == 'amp':
            q = self.session.query(db.PulseResponse.data)
            q = q.join(PulseResponseStrength)
            q = q.filter(PulseResponseStrength.id.in_(ids))
        else:
            q = self.session.query(db.Baseline.data)
            q = q.join(BaselineResponseStrength)
            q = q.filter(BaselineResponseStrength.id.in_(ids))
        recs = q.all()
        
        self.trace_plot.clear()
        self.dec_plot.clear()
        for rec in recs:
            trace = Trace(rec[0], sample_rate=20e3)
            self.trace_plot.plot(trace.time_values, trace.data - np.median(trace.time_slice(1e-3, 9e-3).data))
            dec_trace = deconv_filter(trace)
            self.dec_plot.plot(dec_trace.time_values, dec_trace.data)
        
    def load_conn(self, pair):
        amp_recs = get_amps(self.session, pair)
        base_recs = get_baseline_amps(self.session, pair, limit=len(amp_recs))
        
        by = [rec['pos_dec_amp'] for rec in amp_recs]
        bx = 1 + np.random.random(size=len(by)) * 0.5
        self.amp_sp.setData(bx, by, data=amp_recs, brush=(255, 255, 255, 80))

        ay = [rec['pos_dec_amp'] for rec in base_recs[:len(by)]]
        ax = np.random.random(size=len(ay)) * 0.5
        self.base_sp.setData(ax, ay, data=base_recs[:len(by)], brush=(255, 255, 255, 80))

        dy = [rec['neg_dec_amp'] for rec in amp_recs]
        dx = 3 + np.random.random(size=len(dy)) * 0.5
        self.amp_sp.addPoints(dx, dy, data=amp_recs, brush=(255, 255, 255, 80))
        
        cy = [rec['neg_dec_amp'] for rec in base_recs[:len(dy)]]
        cx = 2 + np.random.random(size=len(cy)) * 0.5
        self.base_sp.addPoints(cx, cy, data=base_recs[:len(dy)], brush=(255, 255, 255, 80))



if __name__ == '__main__':
    #tt = pg.debug.ThreadTrace()


    from multipatch_analysis.ui.multipatch_nwb_viewer import MultipatchNwbViewer    
    from multipatch_analysis.experiment_list import cached_experiments
    expts = cached_experiments()

    pg.dbg()

    if '--rebuild' in sys.argv:
        connection_strength_tables.drop_tables()
        pulse_response_strength_tables.drop_tables()
        init_tables()
        rebuild_strength()
        rebuild_connectivity()
    elif '--rebuild-connectivity' in sys.argv:
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
    
    gl = pg.GraphicsLayoutWidget()
    win.addWidget(gl)
    
    rs_plots = ResponseStrengthPlots(session)
    gl.addItem(rs_plots.str_plot)
    gl.addItem(rs_plots.trace_plot, row=1, col=0)
    gl.addItem(rs_plots.dec_plot, row=2, col=0)
    
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
        sec = time.mktime(ts.timetuple()) + ts.microsecond * 1e-6
        print(sec)

    b.itemSelectionChanged.connect(selected)

    nwb_viewer = MultipatchNwbViewer()

    def dbl_clicked(index):
        item = b.itemFromIndex(index)[0]
        nwb = item.expt.ephys_file
        # temporary workaround:
        nwb = os.path.join(config.cache_path, nwb.split('cache/')[1])
        print(nwb)
        
        nwb_viewer.load_nwb(nwb)
        nwb_viewer.show()
        
    b.doubleClicked.connect(dbl_clicked)

    spw = pg.ScatterPlotWidget()
    spw.style['symbolPen'] = None
    
    spw.show()

    import pandas
    query = """
    select connection_strength.*, ((DATE_PART('day', experiment.acq_timestamp - '1970-01-01'::timestamp) * 24 + 
                DATE_PART('hour', experiment.acq_timestamp - '1970-01-01'::timestamp)) * 60 +
                DATE_PART('minute', experiment.acq_timestamp - '1970-01-01'::timestamp)) * 60 +
                DATE_PART('second', experiment.acq_timestamp - '1970-01-01'::timestamp) as acq_timestamp,
                ABS(connection_strength.amp_med) as abs_amp_med,
                ABS(connection_strength.base_amp_med) as abs_base_amp_med,
                ABS(connection_strength.amp_med_minus_base) as abs_amp_med_minus_base,
                ABS(connection_strength.deconv_amp_med) as abs_deconv_amp_med,
                ABS(connection_strength.deconv_base_amp_med) as abs_deconv_base_amp_med,
                ABS(connection_strength.deconv_amp_med_minus_base) as abs_deconv_amp_med_minus_base,
                experiment.rig_name,
                experiment.acsf,
                pre_cell.cre_type as pre_cre_type,
                pre_cell.target_layer as pre_target_layer,
                post_cell.cre_type as post_cre_type,
                post_cell.target_layer as post_target_layer,
                pair.synapse,
                pair.crosstalk_artifact
    from connection_strength
    join pair on connection_strength.pair_id=pair.id
    join cell pre_cell on pair.pre_cell_id=pre_cell.id
    join cell post_cell on pair.post_cell_id=post_cell.id
    join experiment on pair.expt_id=experiment.id
    """
    df = pandas.read_sql(query, session.bind)
    recs = df.to_records()


    spw.setFields([
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
        ('amp_med', {'units': 'V'}),
        ('abs_amp_med', {'units': 'V'}),
        ('amp_stdev', {'units': 'V'}),
        ('base_amp_med', {'units': 'V'}),
        ('abs_base_amp_med', {'units': 'V'}),
        ('base_amp_stdev', {'units': 'V'}),
        ('amp_med_minus_base', {'units': 'V'}),
        ('abs_amp_med_minus_base', {'units': 'V'}),
        ('amp_stdev_minus_base', {'units': 'V'}),
        ('deconv_amp_med', {'units': 'V'}),
        ('abs_deconv_amp_med', {'units': 'V'}),
        ('deconv_amp_stdev', {'units': 'V'}),
        ('deconv_base_amp_med', {'units': 'V'}),
        ('abs_deconv_base_amp_med', {'units': 'V'}),
        ('deconv_base_amp_stdev', {'units': 'V'}),
        ('deconv_amp_med_minus_base', {'units': 'V'}),
        ('abs_deconv_amp_med_minus_base', {'units': 'V'}),
        ('deconv_amp_stdev_minus_base', {'units': 'V'}),
        ('amp_ttest', {}),
        ('deconv_amp_ttest', {}),
        ('amp_ks2samp', {}),
        ('deconv_amp_ks2samp', {}),
        ('crosstalk_artifact', {'units': 'V'}),
        ('latency_med', {'units': 's'}),
        ('latency_stdev', {'units': 's'}),
        ('base_latency_med', {'units': 's'}),
        ('base_latency_stdev', {'units': 's'}),
    ])

    spw.setData(recs)

    def conn_clicked(spw, points):
        d = points[0].data()
        b.select(pair_id=d['pair_id'])

    spw.sigScatterPlotClicked.connect(conn_clicked)
