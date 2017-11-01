"""
Big question: what's the best way to measure synaptic strength / connectivity?

"""
from __future__ import print_function, division

from collections import OrderedDict
from constants import EXCITATORY_CRE_TYPES, INHIBITORY_CRE_TYPES

import argparse
import sys
import pyqtgraph as pg
import os
import pickle
import io
import multiprocessing
import numpy as np
import scipy.stats

from sqlalchemy.orm import aliased

from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.filter import bessel_filter
from neuroanalysis.event_detection import exp_deconvolve

import database as db



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
    schemas = {
        'pulse_response_strength': [
            ('pulse_response_id', 'pulse_response.id', '', {'index': True}),
            ('pos_amp', 'float'),
            ('neg_amp', 'float'),
            ('pos_base_amp', 'float'),
            ('neg_base_amp', 'float'),
            ('pos_dec_amp', 'float'),
            ('neg_dec_amp', 'float'),
            ('pos_dec_base_amp', 'float'),
            ('neg_dec_base_amp', 'float'),
        ]
        #'deconv_pulse_response': [
            #"Exponentially deconvolved pulse responses",
        #],
    }


    def create_mappings(self):
        TableGroup.create_mappings(self)
        
        PulseResponseStrength = self['pulse_response_strength']
        
        db.PulseResponse.pulse_response_strength = db.relationship(PulseResponseStrength, back_populates="pulse_response", cascade="delete", single_parent=True)
        PulseResponseStrength.pulse_response = db.relationship(db.PulseResponse, back_populates="pulse_response_strength")


class ConnectionStrengthTableGroup(TableGroup):
    schemas = {
        'connection_strength': [
            ('experiment_id', 'experiment.id', '', {'index': True}),
            ('pre_id', 'int', '', {'index': True}),
            ('post_id', 'int', '', {'index': True}),
            ('user_connected', 'bool', 'Whether the experimenter marked this pair as connected.'),
            ('synapse_type', 'str', '"ex" or "in"'),
            ('n_samples', 'int'),
            ('amp_mean', 'float'),
            ('amp_stdev', 'float'),
            ('amp_vom', 'float'),
            ('base_amp_mean', 'float'),
            ('base_amp_stdev', 'float'),
            ('base_amp_vom', 'float'),
            ('deconv_amp_mean', 'float'),
            ('deconv_amp_stdev', 'float'),
            ('deconv_amp_vom', 'float'),
            ('deconv_base_amp_mean', 'float'),
            ('deconv_base_amp_stdev', 'float'),
            ('deconv_base_amp_vom', 'float'),
            ('amp_ttest', 'float'),
            ('deconv_amp_ttest', 'float'),
            ('amp_ks2samp', 'float'),
            ('deconv_amp_ks2samp', 'float'),
        ],
    }


pulse_response_strength_tables = PulseResponseStrengthTableGroup()
connection_strength_tables = ConnectionStrengthTableGroup()

def init_tables():
    global PulseResponseStrength, ConnectionStrength
    pulse_response_strength_tables.create_tables()
    connection_strength_tables.create_tables()

    PulseResponseStrength = pulse_response_strength_tables['pulse_response_strength']
    ConnectionStrength = connection_strength_tables['connection_strength']

def measure_peak(trace, sign, baseline=(0e-3, 9e-3), response=(11e-3, 17e-3)):
    baseline = trace.time_slice(*baseline).data.mean()
    if sign == '+':
        peak = trace.time_slice(*response).data.max()
    else:
        peak = trace.time_slice(*response).data.min()
    return peak - baseline


def measure_sum(trace, sign, baseline=(0e-3, 9e-3), response=(12e-3, 17e-3)):
    baseline = trace.time_slice(*baseline).data.sum()
    peak = trace.time_slice(*response).data.sum()
    return peak - baseline

        
def deconv_filter(trace, tau=15e-3, lowpass=300.):
    dec = exp_deconvolve(trace, tau)
    baseline = np.median(dec.time_slice(trace.t0, trace.t0+10e-3).data)
    deconv = bessel_filter(dec-baseline, lowpass)
    return deconv


@db.default_session
def rebuild_strength(parallel=True, workers=6, session=None):
    print("Rebuilding response strength table..")
    
    max_pulse_id = session.execute('select max(id) from pulse_response').fetchone()[0]
    chunk = 1 + (max_pulse_id // workers)
    parts = [(chunk*i, chunk*(i+1)) for i in range(4)]

    if parallel:
        pool = multiprocessing.Pool(processes=workers)
        pool.map(compute_strength, parts)
    else:
        for part in parts:
            compute_strength(part)

            
@db.default_session
def compute_strength(inds, session=None):
    start_id, stop_id = inds
    q = """SELECT 
        pulse_response.id AS pulse_response_id,
        pulse_response.data AS pulse_response_data,
        baseline.data AS baseline_data 
    FROM 
        pulse_response 
        JOIN baseline ON baseline.id = pulse_response.baseline_id
    WHERE pulse_response.id >= %d and pulse_response.id < %d
    ORDER BY pulse_response.id
    LIMIT 1000
    """
    
    prof = pg.debug.Profiler(disabled=True, delayed=False)
    
    next_id = start_id
    while True:
        response = ses.execute(q % (next_id, stop_id))  # bottleneck here ~30 ms
        prof('exec')
        recs = response.fetchall()
        prof('fetch')
        if len(recs) == 0:
            break
        new_recs = []
        for rec in recs:
            pr_id, data, baseline = rec
            data = np.load(io.BytesIO(data))
            baseline = np.load(io.BytesIO(baseline))
            #prof('load')
            
            new_rec = {'pulse_response_id': pr_id}
            
            data = Trace(data, sample_rate=20e3)
            new_rec['pos_amp'] = measure_peak(data, '+')
            new_rec['neg_amp'] = measure_peak(data, '-')
            #prof('comp 1')
            
            base = Trace(baseline, sample_rate=20e3)
            new_rec['pos_base_amp'] = measure_peak(base, '+')
            new_rec['neg_base_amp'] = measure_peak(base, '-')
            #prof('comp 2')

            dec_data = deconv_filter(data)
            new_rec['pos_dec_amp'] = measure_peak(dec_data, '+')
            new_rec['neg_dec_amp'] = measure_peak(dec_data, '-')
            #prof('comp 3')
            
            dec_base = deconv_filter(base)
            new_rec['pos_dec_base_amp'] = measure_peak(dec_base, '+')
            new_rec['neg_dec_base_amp'] = measure_peak(dec_base, '-')
            #prof('comp 4')
            new_recs.append(new_rec)
        
        next_id = pr_id + 1
        
        sys.stdout.write("%d / %d\r" % (next_id-start_id, stop_id-start_id))
        sys.stdout.flush()

        prof('process')
        ses.bulk_insert_mappings(PulseResponseStrength, new_recs)
        prof('insert')
        new_recs = []

        ses.commit()
        prof('commit')

    
@db.default_session
def rebuild_connectivity(session):
    print("Rebuilding connectivity table..")
    
    # cheating:
    import experiment_list
    expt_cache = experiment_list.cached_experiments()
    
    for expt in list_experiments():
        try:
            cached_expt = expt_cache[expt.acq_timestamp]
        except:
            cached_expt = None
        
        for devs in get_experiment_pairs(expt):
            amps = get_amps(session, expt, devs)
            
            conn = ConnectionStrength(experiment_id=expt.id, pre_id=devs[0], post_id=devs[1])
            
            # Whether the user marked this as connected
            conn.user_connected = None if cached_expt is None else (devs in cached_expt.connections)
            
            # decide whether to treat this connection as excitatory or inhibitory
            # (probably we can do much better here)
            pos_amp = amps['pos_dec_amp'].mean() - amps['pos_dec_base_amp'].mean()
            neg_amp = amps['neg_dec_amp'].mean() - amps['neg_dec_base_amp'].mean()
            if pos_amp > -neg_amp:
                conn.synapse_type = 'ex'
                pfx = 'pos_'
            else:
                conn.synapse_type = 'in'
                pfx = 'neg_'
            # select out positive or negative amplitude columns
            amp, base_amp = amps[pfx+'amp'], amps[pfx+'base_amp']
            dec_amp, dec_base_amp = amps[pfx+'dec_amp'], amps[pfx+'dec_base_amp']
    
            # compute mean/stdev of samples
            n_samp = len(amps)
            conn.n_samples = n_samp
            if n_samp == 0:
                continue
            conn.amp_mean = amp.mean()
            conn.amp_stdev = amp.std()
            conn.amp_vom = amp.var() / n_samp
            conn.base_amp_mean = base_amp.mean()
            conn.base_amp_stdev = base_amp.std()
            conn.base_amp_vom = base_amp.var() / n_samp
            conn.deconv_amp_mean = dec_amp.mean()
            conn.deconv_amp_stdev = dec_amp.std()
            conn.deconv_amp_vom = dec_amp.var() / n_samp
            conn.deconv_base_amp_mean = dec_base_amp.mean()
            conn.deconv_base_amp_stdev = dec_base_amp.std()
            conn.deconv_base_amp_vom = dec_base_amp.var() / n_samp
            
            # do some statistical tests
            conn.amp_ks2samp = scipy.stats.ks_2samp(amp, base_amp).pvalue
            conn.deconv_amp_ks2samp = scipy.stats.ks_2samp(dec_amp, dec_base_amp).pvalue
            conn.amp_ttest = scipy.stats.ttest_ind(amp, base_amp, equal_var=False).pvalue
            conn.deconv_amp_ttest = scipy.stats.ttest_ind(dec_amp, dec_base_amp, equal_var=False).pvalue
            session.add(conn)
        
        session.commit()

@db.default_session
def list_experiments(session):
    return session.query(db.Experiment).all()


@db.default_session
def get_experiment_devs(expt, session):
    devs = session.query(db.Recording.device_key).join(db.SyncRec).filter(db.SyncRec.experiment_id==expt.id)
    
    # Only keep channels with >2 recordings
    counts = {}
    for r in devs:
        counts.setdefault(r[0], 0)
        counts[r[0]] += 1
    devs = [k for k,v in counts.items() if v > 2]
    devs.sort()
    
    return devs


def get_experiment_pairs(expt):
    devs = get_experiment_devs(expt)
    return [(pre, post) for pre in devs for post in devs if pre != post]


class ExperimentBrowser(pg.TreeWidget):
    def __init__(self):
        pg.TreeWidget.__init__(self)
        self.setColumnCount(6)
        self.populate()
        
    def populate(self):
        # cheating:
        import experiment_list
        expts = experiment_list.cached_experiments()
        
        self.session = db.Session()
        for expt in list_experiments(session=self.session):
            date = expt.acq_timestamp
            date_str = date.strftime('%Y-%m-%d')
            slice = expt.slice
            expt_item = pg.TreeWidgetItem(map(str, [date_str, expt.rig_name, slice.species, expt.target_region, slice.genotype]))
            expt_item.expt = expt
            self.addTopLevelItem(expt_item)

            pairs = get_experiment_pairs(expt)
            
            for d1, d2 in pairs:
                pre_id = d1+1
                post_id = d1+2
                
                try:
                    e = expts[date]
                    conn = (pre_id, post_id) in e.connections
                except:
                    conn = 'no expt'
                
                pair_item = pg.TreeWidgetItem(['%d => %d' % (d1, d2), str(conn)])
                expt_item.addChild(pair_item)
                pair_item.devs = (d1, d2)
                pair_item.expt = expt


def get_amps(session, expt, devs, clamp_mode='ic'):
    """Select records from pulse_response_strength table
    """
    q = session.query(
        PulseResponseStrength.id,
        PulseResponseStrength.pos_amp,
        PulseResponseStrength.neg_amp,
        PulseResponseStrength.pos_base_amp,
        PulseResponseStrength.neg_base_amp,
        PulseResponseStrength.pos_dec_amp,
        PulseResponseStrength.neg_dec_amp,
        PulseResponseStrength.pos_dec_base_amp,
        PulseResponseStrength.neg_dec_base_amp,
    ).join(db.PulseResponse)

    dtype = [
        ('id', 'int'),
        ('pos_amp', 'float'),
        ('neg_amp', 'float'),
        ('pos_base_amp', 'float'),
        ('neg_base_amp', 'float'),
        ('pos_dec_amp', 'float'),
        ('neg_dec_amp', 'float'),
        ('pos_dec_base_amp', 'float'),
        ('neg_dec_base_amp', 'float'),
    ]
    
    q, pre_rec, post_rec = join_pulse_response_to_expt(q)
        
    filters = [
        (db.Experiment.id==expt.id,),
        (pre_rec.device_key==devs[0],),
        (post_rec.device_key==devs[1],),
        (db.PatchClampRecording.clamp_mode==clamp_mode,),
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


if __name__ == '__main__':
    import cProfile
    p = cProfile.Profile()
    p.enable()
    
    try:
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
    finally:
        p.disable()
        #p.print_stats(sort='cumulative')
    
    win = pg.QtGui.QSplitter()
    win.setOrientation(pg.QtCore.Qt.Horizontal)
    win.resize(1000, 800)
    win.show()
    
    b = ExperimentBrowser()
    win.addWidget(b)
    
    gl = pg.GraphicsLayoutWidget()
    win.addWidget(gl)
    
    plt = gl.addPlot()

    sp = pg.ScatterPlotItem(pen=None, symbol='o', symbolPen=None)
    plt.addItem(sp)
    
    plt2 = gl.addPlot(row=1, col=0, labels={'bottom': ('time', 's'), 'left': ('Vm', 'V')})

    session = db.Session()
    
    def selected(*args):
        global session
        try:
            sel = b.selectedItems()[0]
            expt = sel.expt
            devs = sel.devs
            
            prof = pg.debug.Profiler(disabled=False)
            recs = get_amps(session, expt, devs)
            prof('query')
            
            ay = [rec['pos_dec_base_amp'] for rec in recs]
            prof('ay')
            ax = np.random.random(size=len(ay)) * 0.5
            by = [rec['pos_dec_amp'] for rec in recs]
            prof('by')
            bx = 1 + np.random.random(size=len(by)) * 0.5
            cy = [rec['neg_dec_base_amp'] for rec in recs]
            prof('cy')
            cx = 2 + np.random.random(size=len(cy)) * 0.5
            dy = [rec['neg_dec_amp'] for rec in recs]
            prof('dy')
            dx = 3 + np.random.random(size=len(dy)) * 0.5
            sp.setData(ax, ay, data=recs, brush=(255, 255, 255, 80))
            sp.addPoints(bx, by, data=recs, brush=(255, 255, 255, 80))
            sp.addPoints(cx, cy, data=recs, brush=(255, 255, 255, 80))
            sp.addPoints(dx, dy, data=recs, brush=(255, 255, 255, 80))
            prof('plot')
        finally:
            session.close()
        

    b.itemSelectionChanged.connect(selected)
    
    def clicked(sp, points):
        global clicked_points, session
        clicked_points = points
        
        ids = [p.data()['id'] for p in points]
        q = session.query(db.PulseResponse.data, db.Baseline.data)
        q = q.join(db.Baseline)
        q = q.join(PulseResponseStrength)
        q = q.filter(PulseResponseStrength.id.in_(ids))
        recs = q.all()
        
        plt2.clear()
        for rec in recs:
            trace = Trace(rec[0], sample_rate=20e3)
            dec_trace = deconv_filter(trace)
            plt2.plot(dec_trace.time_values, dec_trace.data)
            
            base = Trace(rec[1], sample_rate=20e3)
            dec_base = deconv_filter(base)
            plt2.plot(dec_base.time_values, dec_base.data, pen='r')
        
    sp.sigClicked.connect(clicked)
