# coding: utf8
"""
2018 E-E manuscript fig 3:
Analysis of detection limits vs synaptic strength, kinetics, and background noise
"""
from __future__ import print_function, division
from datetime import datetime
import sys
import multiprocessing
import pyqtgraph as pg
import numpy as np
from scipy import stats

from neuroanalysis.data import TSeries, TSeriesList
import strength_analysis
from strength_analysis import TableGroup
from aisynphys.database import database as db


class DetectionLimitTableGroup(TableGroup):
    """Measures pulse amplitudes for each pulse response and background chunk.
    """
    schemas = {
        'detection_limit': [
            ('pair_id', 'pair.id', '', {'index': True}),
            ('simulation_results', 'object'),
            ('minimum_amplitude', 'float'),
        ]
    }

    def create_mappings(self):
        TableGroup.create_mappings(self)
        
        DetectionLimit = self['detection_limit']
        
        db.Pair.detection_limit = db.relationship(DetectionLimit, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
        DetectionLimit.pair = db.relationship(db.Pair, back_populates="detection_limit", single_parent=True)

detection_limit_tables = DetectionLimitTableGroup()

def init_tables():
    global DetectionLimit
    detection_limit_tables.create_tables()
    DetectionLimit = detection_limit_tables['detection_limit']


def measure_limit(pair, session, classifier):
    pair_id = pair.experiment.acq_timestamp, pair.pre_cell.ext_id, pair.post_cell.ext_id
    print(pair_id)

    amps = strength_analysis.get_amps(session, pair)
    base_amps = strength_analysis.get_baseline_amps(session, pair, amps=amps, clamp_mode='ic')
    
    q = strength_analysis.response_query(session)
    q = q.join(strength_analysis.PulseResponseStrength)
    q = q.filter(strength_analysis.PulseResponseStrength.id.in_(amps['id']))
    q = q.join(db.MultiPatchProbe)
    q = q.filter(db.MultiPatchProbe.induction_frequency < 100)

    fg_recs = q.all()

    traces = []
    deconvs = []
    #for rec in fg_recs[:100]:
        #result = strength_analysis.analyze_response_strength(rec, source='pulse_response', lpf=True, lowpass=2000,
                                            #remove_artifacts=False, bsub=True)
        #trace = result['raw_trace']
        #trace.t0 = -result['spike_time']
        #trace = trace - np.median(trace.time_slice(-0.5e-3, 0.5e-3).data)
        #traces.append(trace)            
        #trace_plot.plot(trace.time_values, trace.data, pen=(0, 0, 0, 20))

        #trace = result['dec_trace']
        #trace.t0 = -result['spike_time']
        #trace = trace - np.median(trace.time_slice(-0.5e-3, 0.5e-3).data)
        #deconvs.append(trace)            
        #deconv_plot.plot(trace.time_values, trace.data, pen=(0, 0, 0, 20))

    ## plot average trace
    #mean = TSeriesList(traces).mean()
    #trace_plot.plot(mean.time_values, mean.data, pen={'color':'g', 'width': 2}, shadowPen={'color':'k', 'width': 3}, antialias=True)
    #mean = TSeriesList(deconvs).mean()
    #deconv_plot.plot(mean.time_values, mean.data, pen={'color':'g', 'width': 2}, shadowPen={'color':'k', 'width': 3}, antialias=True)

    #bins = np.arange(-0.001, 0.015, 0.0005) 
    #field = 'pos_dec_amp'
    #n = min(len(amps), len(base_amps))
    #hist_y, hist_bins = np.histogram(base_amps[:n][field], bins=bins)
    #hist_plot.plot(hist_bins, hist_y, stepMode=True, pen=None, brush=(200, 0, 0, 150), fillLevel=0)
    #hist_y, hist_bins = np.histogram(amps[:n][field], bins=bins)
    #hist_plot.plot(hist_bins, hist_y, stepMode=True, pen='k', brush=(0, 150, 150, 100), fillLevel=0)

    q = strength_analysis.baseline_query(session)
    q = q.join(strength_analysis.BaselineResponseStrength)
    q = q.filter(strength_analysis.BaselineResponseStrength.id.in_(base_amps['id']))
    bg_recs = q.all()

    # measure background connection strength
    bg_results = [strength_analysis.analyze_response_strength(rec, 'baseline') for rec in bg_recs]
    bg_results = strength_analysis.str_analysis_result_table(bg_results, bg_recs)

    # for this example, we use background data to simulate foreground
    # (but this will be biased due to lack of crosstalk in background data)
    fg_recs = bg_recs

    # now measure foreground simulated under different conditions
    amps = 2e-6 * 2**np.arange(9)
    amps[0] = 0
    rtime = 2e-3
    dt = 1 / db.default_sample_rate
    results = np.empty(len(amps), dtype=[('results', object), ('prediction', float), ('confidence', float), ('traces', object), ('rise_time', float), ('amp', float)])
    print("  Simulating synaptic events..")
    # pool = multiprocessing.Pool(4)
    limit_entry = DetectionLimit(pair=pair)

    results = []
    avg_conf = []
    limit = None
    for i,amp in enumerate(amps):
        print("----- %d/%d  %0.3g      \r" % (i, len(amps), amp),)
        result = strength_analysis.simulate_connection(fg_recs, bg_results, classifier, amp, rtime)
        results.append({'amp': amp, 'rise_time': rtime, 'predictions': list(result['predictions']), 'confidence': list(result['confidence'])})
        avg_conf.append(result['confidence'].mean())
        print(results[-1])
        # if we crossed threshold, interpolate to estimate the minimum amplitude
        # and terminate the loop early
        if limit is None and i > 0 and avg_conf[-1] > classifier.prob_threshold:
            a1 = amps[i-1]
            a2 = amp
            c1,c2 = avg_conf[-2:]
            s = (classifier.prob_threshold - c1) / (c2 - c1)
            limit = a1 + s * (a2 - a1)
            break

    limit_entry.simulation_results = results
    limit_entry.minimum_amplitude = limit

    session.add(limit_entry)
    session.commit()


def build_detection_limits():
    # silence warnings about fp issues
    np.seterr(all='ignore')

    # read all pair records from DB
    classifier = strength_analysis.get_pair_classifier(seed=0, use_vc_features=False)
    conns = strength_analysis.query_all_pairs(classifier)

    # filter
    mask = np.isfinite(conns['ic_deconv_amp_mean'])
    filtered = conns[mask]

    # remove recordings with gain errors
    mask = filtered['ic_deconv_amp_mean'] < 0.02

    # remove recordings with high crosstalk
    mask &= abs(filtered['ic_crosstalk_mean']) < 60e-6

    # remove recordings with low sample count
    mask &= filtered['ic_n_samples'] > 50

    typs = filtered['pre_cre_type']
    mask &= typs == filtered['post_cre_type']

    typ_mask = ((typs == 'sim1') | (typs == 'tlx3') | (typs == 'unknown') | (typs == 'rorb') | (typs == 'ntsr1'))
    mask &= typ_mask

    filtered = filtered[mask]

    c_mask = filtered['synapse'] == True
    u_mask = ~c_mask

    signal = filtered['confidence']
    background = filtered['ic_base_deconv_amp_mean']

    session = db.session()

    # do selected connections first
    count = 0
    for i,rec in enumerate(filtered):
        print("================== %d/%d ===================== " % (i, len(filtered)))
        pair = session.query(db.Pair).filter(db.Pair.id==rec['pair_id']).all()[0]
        if pair.detection_limit is not None:
            print("    skip!")
            continue
        try:
            measure_limit(pair, session, classifier)
        except Exception:
            sys.excepthook(*sys.exc_info())
        
        count += 1
        if count > 100:
            print("Bailing out before memory fills up.")
            sys.exit(0)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--rebuild', action='store_true', default=False)
    
    args, extra = parser.parse_known_args(sys.argv[1:])

    pg.dbg()
    if args.rebuild and raw_input("Drop and rebuild detection limit table? ") == 'y':
        detection_limit_tables.drop_tables()
        init_tables()
    else:
        init_tables()

    build_detection_limits()

