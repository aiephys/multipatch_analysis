# coding: utf8
"""
2018 E-E manuscript fig 3:
Analysis of detection limits vs synaptic strength, kinetics, and background noise
"""
from __future__ import print_function, division
from datetime import datetime
import pyqtgraph as pg
import numpy as np

from strength_analysis import query_all_pairs


if __name__ == '__main__':
    show_conns = [
        # (expt_uid, pre_cell, post_cell)
        ("low signal, low noise", (1499277786.89, 1, 3)),
        ("low signal, high noise", (1495833911.11, 1, 8)),
        ("high signal, high noise", (1489441647.6, 8, 5)),
    ]
    
    pg.mkQApp()
    win = pg.GraphicsLayoutWidget()
    win.show()
    win.resize(1600, 600)

    scatter_plot = win.addPlot(0, 0, rowspan=len(show_conns))
    scatter_plot.setLogMode(True, True)
    scatter_plot.setAspectLocked()
    scatter_plot.setFixedWidth(600)

    # read all pair records from DB
    conns = query_all_pairs()

    # filter
    mask = np.isfinite(conns['abs_deconv_base_amp_med'])
    filtered = conns[mask]
    # remove recordings with gain errors
    mask = filtered['abs_deconv_base_amp_med'] < 0.02
    # remove recordings likely to have high crosstalk
    mask &= filtered['electrode_distance'] > 1
    # remove recordings with low sample count
    mask &= filtered['n_samples'] > 50

    # filtered = filtered[mask]


    # plot signal vs background for all pairs
    brushes = [pg.mkBrush('y') if c['synapse'] else pg.mkBrush(0.5) for c in filtered]

    c_mask = filtered['synapse']
    u_mask = ~c_mask

    signal = filtered['abs_deconv_amp_med']
    background = filtered['abs_deconv_base_amp_med']

    c_plt = scatter_plot.plot(background[c_mask], signal[c_mask], pen=None, symbol='d', symbolPen=None, symbolBrush='y', symbolSize=8)
    u_plt = scatter_plot.plot(background[u_mask], signal[u_mask], pen=None, symbol='o', symbolPen=None, symbolBrush=(200, 200, 200, 100), symbolSize=4)
    u_plt.setZValue(-10)
    # u_plt.scatter.setCompositionMode(pg.QtGui.QPainter.CompositionMode_Plus)


    for i, c in enumerate(show_conns):
        name, key = c
        ts, pre_id, post_id = key
        idx = np.argwhere((abs(filtered['acq_timestamp'] - ts) < 1) & (filtered['pre_cell_id'] == pre_id) & (filtered['post_cell_id'] == post_id))
        assert idx.size == 1
        idx = idx[0,0]

        scatter_plot.plot([background[idx]], [signal[idx]], pen=None, symbol='o', size=10)

        trace_plt = win.addPlot(i, 1)
        hist_plt = win.addPlot(i, 2)
        limit_plt = win.addPlot(i, 3)







