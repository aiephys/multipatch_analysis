# coding: utf8
"""
2018 E-E manuscript fig 3:
Analysis of detection limits vs synaptic strength, kinetics, and background noise
"""
from __future__ import print_function, division
from datetime import datetime
import pyqtgraph as pg
import numpy as np
from scipy import stats

from neuroanalysis.data import Trace, TraceList
from neuroanalysis.fitting import Psp
import strength_analysis
from multipatch_analysis.database import database as db


if __name__ == '__main__':
    # silence warnings about fp issues
    np.seterr(all='ignore')

    show_conns = [
        # (expt_uid, pre_cell, post_cell)
        # ("low signal, low noise", (1499277786.89, 1, 3)),
        # ("low signal, low noise", (1503953399.15, 1, 7)),
        # ("low signal, high noise", (1495833911.11, 1, 8)),
        # ("low signal, high noise", (1509566559.2, 7, 1)),
        # ("high signal, high noise", (1489441647.6, 8, 5)),

        # ("low signal, low noise", (1506537287.63, 7, 8)),
        ("high signal, low noise", (1499725138.07, 7, 4)),
        # ("low signal, high noise", (1494969844.93, 6, 1)),
        # ("low signal, low noise", (1502312765.01, 1, 4)),
        # ("high noise, no connection", (1489009391.46, 3, 5)),
    ]
    
    pg.mkQApp()
    pg.dbg()
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')

    win = pg.GraphicsLayoutWidget()
    win.show()
    win.resize(1600, 600)

    scatter_plot = win.addPlot(0, 0, rowspan=len(show_conns))
    scatter_plot.setLogMode(True, True)
    scatter_plot.setAspectLocked()
    scatter_plot.setFixedWidth(500)

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


    # plot signal vs background for all pairs
    brushes = [pg.mkBrush('y') if c['synapse'] else pg.mkBrush(0.5) for c in filtered]

    c_mask = filtered['synapse'] == True
    u_mask = ~c_mask

    signal = filtered['confidence']
    background = filtered['ic_base_deconv_amp_mean']

    c_plot = scatter_plot.plot(background[c_mask], signal[c_mask], pen=None, symbol='d', symbolPen='k', symbolBrush=(0, 255, 255), symbolSize=10, data=filtered[c_mask])

    u_plot = scatter_plot.plot(background[u_mask], signal[u_mask], pen=None, symbol='o', symbolPen=None, symbolBrush=(50, 50, 50, 80), symbolSize=4, data=filtered[u_mask])
    u_plot.setZValue(-10)
    # u_plot.scatter.setCompositionMode(pg.QtGui.QPainter.CompositionMode_Plus)

    trace_plots = []
    deconv_plots = []
    hist_plots = []

    session = db.Session()

    def add_connection_plots(i, name, timestamp, pre_id, post_id):
        global session, win, filtered
        p = pg.debug.Profiler(disabled=True, delayed=False)
        trace_plot = win.addPlot(i, 1)
        trace_plots.append(trace_plot)
        deconv_plot = win.addPlot(i, 2)
        deconv_plots.append(deconv_plot)
        hist_plot = win.addPlot(i, 3)
        hist_plots.append(hist_plot)
        limit_plot = win.addPlot(i, 4)
        limit_plot.addLegend()
        limit_plot.setLogMode(True, False)
        # Find this connection in the pair list
        idx = np.argwhere((abs(filtered['acq_timestamp'] - timestamp) < 1) & (filtered['pre_cell_id'] == pre_id) & (filtered['post_cell_id'] == post_id))
        if idx.size == 0:
            print("not in filtered connections")
            return
        idx = idx[0,0]
        p()

        # Mark the point in scatter plot
        scatter_plot.plot([background[idx]], [signal[idx]], pen='k', symbol='o', size=10, symbolBrush='r', symbolPen=None)
        
        # Plot example traces and histograms
        for plts in [trace_plots, deconv_plots]:
            plt = plts[-1]
            plt.setXLink(plts[0])
            plt.setYLink(plts[0])
            plt.setXRange(-10e-3, 17e-3, padding=0)
            plt.hideAxis('left')
            plt.hideAxis('bottom')
            plt.addLine(x=0)
            plt.setDownsampling(auto=True, mode='peak')
            plt.setClipToView(True)
            hbar = pg.QtGui.QGraphicsLineItem(0, 0, 2e-3, 0)
            hbar.setPen(pg.mkPen(color='k', width=5))
            plt.addItem(hbar)
            vbar = pg.QtGui.QGraphicsLineItem(0, 0, 0, 100e-6)
            vbar.setPen(pg.mkPen(color='k', width=5))
            plt.addItem(vbar)


        hist_plot.setXLink(hist_plots[0])
        
        pair = session.query(db.Pair).filter(db.Pair.id==filtered[idx]['pair_id']).all()[0]
        p()
        amps = strength_analysis.get_amps(session, pair)
        p()
        base_amps = strength_analysis.get_baseline_amps(session, pair)
        p()
        
        q = strength_analysis.response_query(session)
        p()
        q = q.join(strength_analysis.PulseResponseStrength)
        q = q.filter(strength_analysis.PulseResponseStrength.id.in_(amps['id']))
        q = q.join(db.MultiPatchProbe)
        q = q.filter(db.MultiPatchProbe.induction_frequency < 100)
        # pre_cell = db.aliased(db.Cell)
        # post_cell = db.aliased(db.Cell)
        # q = q.join(db.Pair).join(db.Experiment).join(pre_cell, db.Pair.pre_cell_id==pre_cell.id).join(post_cell, db.Pair.post_cell_id==post_cell.id)
        # q = q.filter(db.Experiment.id==filtered[idx]['experiment_id'])
        # q = q.filter(pre_cell.ext_id==pre_id)
        # q = q.filter(post_cell.ext_id==post_id)

        fg_recs = q.all()
        p()

        traces = []
        deconvs = []
        for rec in fg_recs[:100]:
            result = strength_analysis.analyze_response_strength(rec, source='pulse_response', lpf=True, lowpass=2000,
                                                remove_artifacts=False, bsub=True)
            trace = result['raw_trace']
            trace.t0 = -result['spike_time']
            trace = trace - np.median(trace.time_slice(-0.5e-3, 0.5e-3).data)
            traces.append(trace)            
            trace_plot.plot(trace.time_values, trace.data, pen=(0, 0, 0, 20))

            trace = result['dec_trace']
            trace.t0 = -result['spike_time']
            trace = trace - np.median(trace.time_slice(-0.5e-3, 0.5e-3).data)
            deconvs.append(trace)            
            deconv_plot.plot(trace.time_values, trace.data, pen=(0, 0, 0, 20))

        # plot average trace
        mean = TraceList(traces).mean()
        trace_plot.plot(mean.time_values, mean.data, pen={'color':'g', 'width': 2}, shadowPen={'color':'k', 'width': 3}, antialias=True)
        mean = TraceList(deconvs).mean()
        deconv_plot.plot(mean.time_values, mean.data, pen={'color':'g', 'width': 2}, shadowPen={'color':'k', 'width': 3}, antialias=True)

        # add label
        label = pg.LabelItem(name)
        label.setParentItem(trace_plot)


        p("analyze_response_strength")

        # bins = np.arange(-0.0005, 0.002, 0.0001) 
        # field = 'pos_amp'
        bins = np.arange(-0.001, 0.015, 0.0005) 
        field = 'pos_dec_amp'
        n = min(len(amps), len(base_amps))
        hist_y, hist_bins = np.histogram(base_amps[:n][field], bins=bins)
        hist_plot.plot(hist_bins, hist_y, stepMode=True, pen=None, brush=(200, 0, 0, 150), fillLevel=0)
        hist_y, hist_bins = np.histogram(amps[:n][field], bins=bins)
        hist_plot.plot(hist_bins, hist_y, stepMode=True, pen='k', brush=(0, 150, 150, 100), fillLevel=0)
        p()

        pg.QtGui.QApplication.processEvents()


        # Plot detectability analysis
        q = strength_analysis.baseline_query(session)
        q = q.join(strength_analysis.BaselineResponseStrength)
        q = q.filter(strength_analysis.BaselineResponseStrength.id.in_(base_amps['id']))
        # q = q.limit(100)
        bg_recs = q.all()

        def clicked(sp, pts):
            data = pts[0].data()
            print("-----------------------\nclicked:", data['rise_time'], data['amp'], data['prediction'], data['confidence'])
            for r in data['results']:
                print({k:r[k] for k in classifier.features})
            traces = data['traces']
            plt = pg.plot()
            bsub = [t.copy(data=t.data - np.median(t.time_slice(0, 1e-3).data)) for t in traces]
            for t in bsub:
                plt.plot(t.time_values, t.data, pen=(0, 0, 0, 50))
            mean = TraceList(bsub).mean()
            plt.plot(mean.time_values, mean.data, pen='g')


        def analyze_response_strength(recs, source, dtype):
            """Wraps strength_analysis.analyze_response_strength to look like
            the result was queried from the DB using get_amps() or get_baseline()
            """
            results = np.empty(len(recs), dtype=dtype)
            for i,rec in enumerate(recs):
                result = strength_analysis.analyze_response_strength(rec, source)
                for key in ['ex_qc_pass', 'in_qc_pass', 'clamp_mode']:
                    result[key] = getattr(rec, key)
                for key,val in result.items():
                    if key in results.dtype.names:
                        results[i][key] = val
            return results

        # measure background connection strength
        bg_results = analyze_response_strength(bg_recs, 'baseline', base_amps.dtype)

        # for this example, we use background data to simulate foreground
        # (but this will be biased due to lack of crosstalk in background data)
        fg_recs = bg_recs

        # now measure foreground simulated under different conditions
        amps = 5e-6 * 1.55**np.arange(6)
        amps[0] = 0
        rtimes = [1e-3, 2e-3, 4e-3, 6e-3]
        dt = 1 / db.default_sample_rate
        results = np.empty((len(amps), len(rtimes)), dtype=[('results', object), ('prediction', bool), ('confidence', float), ('traces', object), ('rise_time', float), ('amp', float)])
        print("  Simulating synaptic events..")
        for j,rtime in enumerate(rtimes):
            for i,amp in enumerate(amps):
                print("---------------------------------------    %d/%d  %d/%d      \r" % (i,len(amps),j,len(rtimes)),)
                t = np.arange(0, 15e-3, dt)
                template = Psp.psp_func(t, xoffset=0, yoffset=0, rise_time=rtime, decay_tau=15e-3, amp=1, rise_power=2)

                results[i,j]['results'] = []
                results[i,j]['rise_time'] = rtime
                results[i,j]['amp'] = amp
                for k in range(6):
                    r_amps = stats.binom.rvs(p=0.2, n=24, size=len(fg_recs)) * stats.norm.rvs(scale=0.3, loc=1, size=len(fg_recs))
                    r_amps *= amp / r_amps.mean()
                    r_latency = np.random.normal(size=len(fg_recs), scale=200e-6, loc=13e-3)
                    fg_results = np.empty(len(fg_recs), dtype=bg_results.dtype)
                    traces = []
                    for k,rec in enumerate(fg_recs):
                        data = rec.data.copy()
                        start = int(r_latency[k] / dt)
                        length = len(rec.data) - start
                        rec.data[start:] += template[:length] * r_amps[k]

                        fg_result = strength_analysis.analyze_response_strength(rec, 'baseline')
                        for key in ['ex_qc_pass', 'in_qc_pass', 'clamp_mode']:
                            fg_result[key] = getattr(rec, key)
                        for key,val in fg_result.items():
                            if key in fg_results.dtype.names:
                                fg_results[k][key] = val

                        traces.append(Trace(rec.data.copy(), dt=dt))
                        traces[-1].amp = r_amps[k]
                        rec.data[:] = data  # can't modify rec, so we have to muck with the array (and clean up afterward) instead
                        
                    conn_result = strength_analysis.analyze_pair_connectivity({('ic', 'fg'): fg_results, ('ic', 'bg'): bg_results, ('vc', 'fg'): [], ('vc', 'bg'): []}, sign=1)
                    results[i,j]['results'].append(conn_result)
                    results[i,j]['traces'] = traces[:100]
                    print(".",)
                    # print(conn_result)
                    print(dict([(k, conn_result[k]) for k in classifier.features]))

                pred = classifier.predict(results[i,j]['results'])
                results[i,j]['prediction'] = pred['prediction'].mean()
                results[i,j]['confidence'] = pred['confidence'].mean()
                print("\nrise time:", rtime, " amplitude:", amp)
                print(pred)


            # c = limit_plot.plot(rtimes, results[i]['result'], pen=(i, len(amps)*1.3), symbol='o', antialias=True, name="%duV"%(amp*1e6), data=results[i], symbolSize=4)
            # c.scatter.sigClicked.connect(clicked)
            # pg.QtGui.QApplication.processEvents()
            c = limit_plot.plot(amps, results[:,j]['confidence'], pen=(j, len(rtimes)*1.3), symbol='o', antialias=True, name="%dus"%(rtime*1e6), data=results[:,j], symbolSize=4)
            c.scatter.sigClicked.connect(clicked)
            pg.QtGui.QApplication.processEvents()

                
        pg.QtGui.QApplication.processEvents()

    # Handle selected individual connections
    next_row = 0
    for name, key in show_conns:
        ts, pre_id, post_id = key
        add_connection_plots(next_row, name, ts, pre_id, post_id)
        next_row += 1

    def clicked(sp, pts):
        global next_row
        print("-----------")
        for pt in pts:
            d = pt.data()
            print(d['acq_timestamp'], d['pre_cell_id'], d['post_cell_id'])
            add_connection_plots(next_row, "", d['acq_timestamp'], d['pre_cell_id'], d['post_cell_id'])
            next_row += 1

    c_plot.scatter.sigClicked.connect(clicked)
    u_plot.scatter.sigClicked.connect(clicked)
