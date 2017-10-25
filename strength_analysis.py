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
    schemas = {
        'pulse_response_strength': [
            ('pulse_response_id', 'pulse_response.id'),
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

    def __init__(self):
        self.mappings = {}
        self.create_mappings()

    def __getitem__(self, item):
        return self.mappings[item]

    def create_mappings(self):
        for k,schema in self.schemas.items():
            self.mappings[k] = db.generate_mapping(k, schema)

        PulseResponseStrength = self['pulse_response_strength']
        
        db.PulseResponse.pulse_response_strength = db.relationship(PulseResponseStrength, back_populates="pulse_response", cascade="delete", single_parent=True)
        PulseResponseStrength.pulse_response = db.relationship(db.PulseResponse, back_populates="pulse_response_strength")

    def drop_tables(self):
        for k in self.schemas:
            if k in db.engine.table_names():
                self[k].__table__.drop(bind=db.engine)

    def create_tables(self):
        for k in self.schemas:
            if k not in db.engine.table_names():
                self[k].__table__.create(bind=db.engine)


tables = TableGroup()

if '--reset-db' in sys.argv:
    tables.drop_tables()

tables.create_tables()



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


def rebuild_tables(parallel=True, workers=6):
    tables.drop_tables()
    tables.create_tables()
    
    ses = db.Session()
    max_pulse_id = ses.execute('select max(id) from pulse_response').fetchone()[0]
    chunk = 1 + (max_pulse_id // workers)
    parts = [(chunk*i, chunk*(i+1)) for i in range(4)]

    if parallel:
        pool = multiprocessing.Pool(processes=workers)
        pool.map(compute_strength, parts)
    else:
        for part in parts:
            compute_strength(part)


def compute_strength(inds):
    start_id, stop_id = inds
    ses = db.Session()
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
        ses.bulk_insert_mappings(tables['pulse_response_strength'], new_recs)
        prof('insert')
        new_recs = []

        ses.commit()
        prof('commit')
    


def _rebuild_tables():
    # ORM version (too slow)
    tables.drop_tables()
    tables.create_tables()
    
    ses = db.Session()
    q = ses.query(db.PulseResponse.id, db.PulseResponse.data, db.Baseline.data).join(db.Baseline).yield_per(100)
    for i,rec in enumerate(q):
        pr_id, data, baseline = rec
        data = Trace(data, sample_rate=20e3)
        prs = tables['pulse_response_strength'](pulse_response_id=pr_id)
        prs.pos_amp = measure_peak(data, '+')
        prs.neg_amp = measure_peak(data, '-')
        
        base = Trace(baseline, sample_rate=20e3)
        prs.pos_base_amp = measure_peak(base, '+')
        prs.neg_base_amp = measure_peak(base, '-')

        dec_data = deconv_filter(data)
        prs.pos_dec_amp = measure_peak(dec_data, '+')
        prs.neg_dec_amp = measure_peak(dec_data, '-')
        
        dec_base = deconv_filter(base)
        prs.pos_dec_base_amp = measure_peak(dec_base, '+')
        prs.neg_dec_base_amp = measure_peak(dec_base, '-')
        
        ses.add(prs)
        
        if i%100 == 0:
            print("%d/%d" % (i, q.count()))
            
        ses.commit()


class ExperimentBrowser(pg.TreeWidget):
    def __init__(self):
        pg.TreeWidget.__init__(self)
        self.setColumnCount(4)
        self.populate()
        
    def populate(self):
        s = db.Session()
        for expt in s.query(db.Experiment):
            date = expt.acq_timestamp.strftime('%Y-%m-%d')
            slice = expt.slice
            expt_item = pg.TreeWidgetItem(map(str, [date, slice.species, expt.target_region, slice.genotype]))
            expt_item.expt = expt
            self.addTopLevelItem(expt_item)

            devs = s.query(db.Recording.device_key).join(db.SyncRec).filter(db.SyncRec.experiment_id==expt.id)
            
            # Only keep channels with >2 recordings
            counts = {}
            for r in devs:
                counts.setdefault(r[0], 0)
                counts[r[0]] += 1
            devs = [k for k,v in counts.items() if v > 2]
            devs.sort()
            
            for d1 in devs:
                for d2 in devs:
                    if d1 == d2:
                        continue
                    pair_item = pg.TreeWidgetItem(['%d => %d' % (d1, d2)])
                    expt_item.addChild(pair_item)
                    pair_item.devs = (d1, d2)
                    pair_item.expt = expt


def get_amps(session, expt, devs, clamp_mode='ic'):
    """Select records from pulse_response_strength table
    """
    pre_rec = aliased(db.Recording)
    post_rec = aliased(db.Recording)
    q = session.query(
        tables['pulse_response_strength'].id,
        tables['pulse_response_strength'].pos_amp,
        tables['pulse_response_strength'].neg_amp,
        tables['pulse_response_strength'].pos_base_amp,
        tables['pulse_response_strength'].neg_base_amp,
        tables['pulse_response_strength'].pos_dec_amp,
        tables['pulse_response_strength'].neg_dec_amp,
        tables['pulse_response_strength'].pos_dec_base_amp,
        tables['pulse_response_strength'].neg_dec_base_amp,
    )
    
    joins = [
        (db.PulseResponse,),
        (post_rec, db.PulseResponse.recording),
        (db.PatchClampRecording,),
        (db.SyncRec,),
        (db.Experiment,),
        (db.StimPulse, db.PulseResponse.stim_pulse),
        (pre_rec, db.StimPulse.recording),
    ]
    for join_args in joins:
        q = q.join(*join_args)
        
    filters = [
        (db.Experiment.id==expt.id,),
        (pre_rec.device_key==devs[0],),
        (post_rec.device_key==devs[1],),
        (db.PatchClampRecording.clamp_mode==clamp_mode,),
    ]
    for filter_args in filters:
        q = q.filter(*filter_args)
    
    return q


if __name__ == '__main__':
    pg.dbg()
    if '--rebuild' in sys.argv:
        rebuild_tables()
    
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
        sel = b.selectedItems()[0]
        expt = sel.expt
        devs = sel.devs
        
        prof = pg.debug.Profiler(disabled=False)
        q = get_amps(session, expt, devs)
        prof('query')
        
        recs = q.all()
        prof('fetch')
        
        ay = [rec[7] for rec in recs]
        prof('ay')
        ax = np.random.random(size=len(ay)) * 0.5
        by = [rec[5] for rec in recs]
        prof('by')
        bx = 1 + np.random.random(size=len(by)) * 0.5
        cy = [rec[8] for rec in recs]
        prof('cy')
        cx = 2 + np.random.random(size=len(cy)) * 0.5
        dy = [rec[6] for rec in recs]
        prof('dy')
        dx = 3 + np.random.random(size=len(dy)) * 0.5
        sp.setData(ax, ay, data=recs, brush=(255, 255, 255, 80))
        sp.addPoints(bx, by, data=recs, brush=(255, 255, 255, 80))
        sp.addPoints(cx, cy, data=recs, brush=(255, 255, 255, 80))
        sp.addPoints(dx, dy, data=recs, brush=(255, 255, 255, 80))
        prof('plot')
        

    b.itemSelectionChanged.connect(selected)
    
    
    def clicked(sp, points):
        global clicked_points
        clicked_points = points
        
        plt2.clear()
        for p in points:
            rec = p.data()
            trace = Trace(rec[0].pulse_response.data, sample_rate=20e3)
            dec_trace = deconv_filter(trace)
            plt2.plot(dec_trace.time_values, dec_trace.data)
            
            base = Trace(rec[0].pulse_response.baseline.data, sample_rate=20e3)
            dec_base = deconv_filter(base)
            plt2.plot(dec_base.time_values, dec_base.data, pen='r')
        
    sp.sigClicked.connect(clicked)

#if __name__ == '__main__':
    #app = pg.mkQApp()
    ##pg.dbg()
    
    #expt_index = sys.argv[1]
    #pre_id, post_id = map(int, sys.argv[2:4])
    
    ## Load experiment index
    #cache_file = 'expts_cache.pkl'
    #expts = ExperimentList(cache=cache_file)

    #expt = expts[expt_index]
    
    #sign = '-' if expt.cells[pre_id].cre_type in INHIBITORY_CRE_TYPES else '+'
    #print("sign:", sign)
    ##analyzer = MultiPatchExperimentAnalyzer(expt.data)
    ##pulses = analyzer.get_evoked_responses(pre_id, post_id, clamp_mode='ic', pulse_ids=[0])
    
    #analyzer = DynamicsAnalyzer(expt, pre_id, post_id, align_to='spike')
    
    ## collect all first pulse responses
    ##responses = analyzer.amp_group
    
    ## collect all events
    #responses = analyzer.all_events

    
    #n_responses = len(responses)
    
    ## do exponential deconvolution on all responses
    #deconv = TraceList()
    #grid1 = PlotGrid()
    #grid1.set_shape(2, 1)
    #grid1[0, 0].setLabels(left=('PSP', 'V'))
    #grid1[1, 0].setLabels(bottom=('time', 's'))
    
    #results = OrderedDict()
    #raw_group = EvokedResponseGroup()
    #deconv_group = EvokedResponseGroup()
    
    #if len(responses) == 0:
        #raise Exception("No data found for this synapse")


    #def input_filter(trace):
        #bsub = trace - np.median(trace.time_slice(0, 10e-3).data)
        #filt = bessel_filter(bsub, 1000.)
        #return filt
    
    #def deconv_filter(trace):
        #dec = exp_deconvolve(trace, 15e-3)
        #baseline = np.median(dec.time_slice(trace.t0, trace.t0+10e-3).data)
        #deconv = bessel_filter(dec-baseline, 300.)
        #return deconv
        

    #order = np.argsort([t.start_time for t in responses.responses])
    #with pg.ProgressDialog("Measuring amplitudes...", maximum=n_responses) as dlg:
        #for i in range(n_responses):
            #r = responses.responses[order[i]]
            #b = responses.baselines[order[i]]            
            #r.t0 = 0
            #b.t0 = 0

            #add_to_avg = True
            #stim_name = r.parent.parent.meta['stim_name']
            #if '200Hz' in stim_name or '100Hz' in stim_name:
                #print("skipped in average:", r.parent.parent)
                #add_to_avg = False
            
            ## lowpass raw data
            #filt = input_filter(r)
            #base_filt = input_filter(b)
            #grid1[0, 0].plot(filt.time_values, filt.data, pen=(255, 255, 255, 100))
            #if add_to_avg:
                #raw_group.add(filt, None)
            
            #results.setdefault('raw_peak', []).append((
                #measure_peak(filt, sign),
                #measure_peak(base_filt, sign)
            #))
            
            ##results.setdefault('raw_sum', []).append((
                ##measure_sum(filt, sign),
                ##measure_sum(base_filt, sign)
            ##))
            
            ## deconvolve
            #deconv = deconv_filter(r)
            #base_deconv = deconv_filter(b)
            #grid1[1, 0].plot(deconv.time_values, deconv.data, pen=(255, 255, 255, 100))
            #if add_to_avg:
                #deconv_group.add(deconv, None)
            
            #results.setdefault('deconv_peak', []).append((
                #measure_peak(deconv, sign),
                #measure_peak(base_deconv, sign)
            #))
            
            ##results.setdefault('deconv_sum', []).append((
                ##measure_sum(deconv, sign),
                ##measure_sum(base_deconv, sign)
            ##))
            
            #dlg += 1
            #if dlg.wasCanceled():
                #break
        
        
    
    #grid1.show()
    #raw_mean = raw_group.mean()
    #grid1[0, 0].plot(raw_mean.time_values, raw_mean.data, pen={'color': 'g', 'width': 2}, shadowPen={'color': 'k', 'width': 3}, antialias=True)
    
    #deconv_mean = deconv_group.mean()
    #grid1[1, 0].plot(deconv_mean.time_values, deconv_mean.data, pen={'color': 'g', 'width': 2}, shadowPen={'color': 'k', 'width': 3}, antialias=True)

    
    #plts = PlotGrid()
    #plts.set_shape(1, len(results))
    
    #for i,k in enumerate(results):
        #amps = np.array([x[0] for x in results[k]])
        #bases = np.array([x[1] for x in results[k]])
        #x = np.linspace(0.0, 1.0, len(amps))
        
        #amps = amps / bases.std()
        #bases = bases / bases.std()
        
        #ks_s, ks_p = scipy.stats.ks_2samp(amps, bases)
        #print("%s%sks_2samp: %g p=%g" % (k, ' '*(30-len(k)), ks_s, ks_p))
        
        #plt = plts[0, i]
        #plt.plot(x, amps, pen=None, symbol='o', symbolBrush=(255, 255, 0, 150), symbolPen=None)
        #plt.plot(x, bases, pen=None, symbol='o', symbolBrush=(255, 0, 0, 150), symbolPen=None)
        #plt.setTitle('%s<br>ks: %0.2g %0.2g' % (k, ks_s, ks_p))
        #if i > 0:
            #plt.hideAxis('left')
        
    #plts.setXLink(plts[0,0])
    #plts.setYLink(plts[0,0])
    #plts[0, 0].setXRange(-0.2, 1.2)
    #plts.show()

    

    
    