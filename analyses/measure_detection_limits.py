# coding: utf8
"""
2018 E-E manuscript fig 3:
Analysis of detection limits vs synaptic strength, kinetics, and background noise
"""
from __future__ import print_function, division
from datetime import datetime
import multiprocessing
import pyqtgraph as pg
import numpy as np
from scipy import stats

from neuroanalysis.data import Trace, TraceList
import strength_analysis
from multipatch_analysis.database import database as db



def measure_limit(timestamp, pre_id, post_id):
    # Find this connection in the pair list
    idx = np.argwhere((abs(filtered['acq_timestamp'] - timestamp) < 1) & (filtered['pre_cell_id'] == pre_id) & (filtered['post_cell_id'] == post_id))
    if idx.size == 0:
        print("not in filtered connections")
        return
    idx = idx[0,0]
    
    pair = session.query(db.Pair).filter(db.Pair.id==filtered[idx]['pair_id']).all()[0]
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
    #mean = TraceList(traces).mean()
    #trace_plot.plot(mean.time_values, mean.data, pen={'color':'g', 'width': 2}, shadowPen={'color':'k', 'width': 3}, antialias=True)
    #mean = TraceList(deconvs).mean()
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
    for i,amp in enumerate(amps):
        print("---------------------------------------    %d/%d      \r" % (i, len(amps)),)
        result = strength_analysis.simulate_connection(fg_recs, bg_results, classifier, amp, rtime, pair_id=(timestamp, pre_id, post_id))


if __name__ == '__main__':
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

    session = db.Session()

    # do selected connections first
    for i,rec in enumerate(filtered):
        pair_id = (rec['acq_timestamp'], rec['pre_cell_id'], rec['post_cell_id'])
        print("================== %s %s %s ===================== (%d/%d)" % (pair_id + (i, len(filtered))))
        try:
            measure_limit(*pair_id)
        except Exception:
            sys.excepthook(*sys.exc_info())
            
    # then everything else
    for i,rec in enumerate(conns):
        pair_id = (rec['acq_timestamp'], rec['pre_cell_id'], rec['post_cell_id'])
        print("================== %s %s %s ===================== (%d/%d)" % (pair_id + (i, len(conns)))
        try:
            measure_limit(*pair_id)
        except Exception:
            sys.excepthook(*sys.exc_info())
            
    
