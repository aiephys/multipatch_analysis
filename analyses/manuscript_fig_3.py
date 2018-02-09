# coding: utf8
"""
2018 E-E manuscript fig 3:
Analysis of detection limits vs synaptic strength, kinetics, and background noise
"""
from __future__ import print_function, division
from datetime import datetime
import pyqtgraph as pg
import numpy as np

import strength_analysis
from multipatch_analysis.database import database as db


if __name__ == '__main__':
    show_conns = [
        # (expt_uid, pre_cell, post_cell)
        ("low signal, low noise", (1499277786.89, 1, 3)),
        ("low signal, high noise", (1495833911.11, 1, 8)),
        ("high signal, high noise", (1489441647.6, 8, 5)),
    ]
    
    pg.mkQApp()
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')

    win = pg.GraphicsLayoutWidget()
    win.show()
    win.resize(1600, 600)

    scatter_plot = win.addPlot(0, 0, rowspan=len(show_conns))
    scatter_plot.setLogMode(True, True)
    scatter_plot.setAspectLocked()
    scatter_plot.setFixedWidth(600)

    # read all pair records from DB
    conns = strength_analysis.query_all_pairs()

    # filter
    mask = np.isfinite(conns['abs_deconv_base_amp_med'])
    filtered = conns[mask]
    # remove recordings with gain errors
    mask = filtered['abs_deconv_base_amp_med'] < 0.02
    # remove recordings likely to have high crosstalk
    mask &= filtered['electrode_distance'] > 1
    # remove recordings with low sample count
    mask &= filtered['n_samples'] > 50

    filtered = filtered[mask]


    # plot signal vs background for all pairs
    brushes = [pg.mkBrush('y') if c['synapse'] else pg.mkBrush(0.5) for c in filtered]

    c_mask = filtered['synapse']
    u_mask = ~c_mask

    signal = filtered['abs_deconv_amp_med']
    background = filtered['abs_deconv_base_amp_med']

    def clicked(sp, pts):
        print("-----------")
        for pt in pts:
            d = pt.data()
            print(d['acq_timestamp'], d['pre_cell_id'], d['post_cell_id'])
    c_plot = scatter_plot.plot(background[c_mask], signal[c_mask], pen=None, symbol='d', symbolPen='k', symbolBrush=(0, 255, 255), symbolSize=10, data=filtered[c_mask])
    c_plot.scatter.sigClicked.connect(clicked)

    u_plot = scatter_plot.plot(background[u_mask], signal[u_mask], pen=None, symbol='o', symbolPen=None, symbolBrush=(50, 50, 50, 80), symbolSize=4)
    u_plot.setZValue(-10)
    # u_plot.scatter.setCompositionMode(pg.QtGui.QPainter.CompositionMode_Plus)

    trace_plots = []
    hist_plots = []

    # Handle selected individual connections
    session = db.Session()
    p = pg.debug.Profiler(disabled=False, delayed=False)
    for i, c in enumerate(show_conns):
        name, key = c
        ts, pre_id, post_id = key
        
        # Find this connection in the pair list
        idx = np.argwhere((abs(filtered['acq_timestamp'] - ts) < 1) & (filtered['pre_cell_id'] == pre_id) & (filtered['post_cell_id'] == post_id))
        assert idx.size == 1
        idx = idx[0,0]
        p()

        # Mark the point in scatter plot
        scatter_plot.plot([background[idx]], [signal[idx]], pen='k', symbol='o', size=10, symbolBrush='r', symbolPen=None)
            
        # Plot example traces and histograms
        trace_plot = win.addPlot(i, 1)
        trace_plots.append(trace_plot)
        trace_plot.setXLink(trace_plots[0])
        trace_plot.setYLink(trace_plots[0])
        trace_plot.setXRange(-10e-3, 20e-3)

        hist_plot = win.addPlot(i, 2)
        hist_plots.append(hist_plot)
        hist_plot.setXLink(hist_plots[0])
        
        pair = session.query(db.Pair).filter(db.Pair.id==filtered[idx]['pair_id']).all()[0]
        p()
        amps = strength_analysis.get_amps(session, pair)
        p()
        base_amps = strength_analysis.get_baseline_amps(session, pair)
        p()
        
        ids = amps['id']
        q = strength_analysis.response_query(session)
        p()
        q = q.join(strength_analysis.PulseResponseStrength)
        q = q.filter(strength_analysis.PulseResponseStrength.id.in_(ids))

        # pre_cell = db.aliased(db.Cell)
        # post_cell = db.aliased(db.Cell)
        # q = q.join(db.Pair).join(db.Experiment).join(pre_cell, db.Pair.pre_cell_id==pre_cell.id).join(post_cell, db.Pair.post_cell_id==post_cell.id)
        # q = q.filter(db.Experiment.id==filtered[idx]['experiment_id'])
        # q = q.filter(pre_cell.ext_id==pre_id)
        # q = q.filter(post_cell.ext_id==post_id)

        q = q.limit(50)
        recs = q.all()
        p()
        print(len(recs))
        for rec in recs:
            result = strength_analysis.analyze_response_strength(rec, source='pulse_response', lpf=True, lowpass=2000,
                                                remove_artifacts=False, bsub=True)
            trace = result['raw_trace']
            trace.t0 = -result['spike_time']
            base = np.median(trace.time_slice(-0.5e-3, 0.5e-3).data)
            trace_plot.plot(trace.time_values, trace.data - base, pen=(0, 0, 0, 20))
        
        p("analyze_response_strength")

        bins = np.arange(-0.0005, 0.002, 0.0001) 
        field = 'pos_amp'
        # bins = np.arange(-0.005, 0.02, 0.001) 
        # field = 'pos_dec_amp'
        hist_y, hist_bins = np.histogram(base_amps[field], bins=bins, density=True)
        hist_plot.plot(hist_bins, hist_y, stepMode=True, pen=None, brush=(200, 0, 0, 150), fillLevel=0)
        hist_y, hist_bins = np.histogram(amps[field], bins=bins, density=True)
        hist_plot.plot(hist_bins, hist_y, stepMode=True, pen='k', brush=(0, 150, 150, 150), fillLevel=0)
        p()

        # Plot detectability analysis
        limit_plot = win.addPlot(i, 3)

        # amps = np.arange(50e-6, 500e-6, 50e-6)
        # rtimes = np.arange(1e-3, 2e-3, 0.1e-3)
        # for amp in amps:
        #     for rtime in rtimes:
                


