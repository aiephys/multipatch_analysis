# coding: utf8
"""
2018 E-E manuscript fig 3:
Analysis of detection limits vs synaptic strength, kinetics, and background noise
"""
from __future__ import print_function, division
import pyqtgraph as pg
import numpy as np

from strength_analysis import query_all_connections


if __name__ == '__main__':
    show_conns = [
        # (expt_uid, pre_cell, post_cell)
    ]
    
    pg.mkQApp()
    win = pg.GraphicsLayoutWidget()
    win.show()
    scatter_plot = win.addPlot(0, 0)
    scatter_plot.setLogMode(True, True)
    scatter_plot.setAspectLocked()

    conns = query_all_connections()


    mask = np.isfinite(conns['abs_deconv_base_amp_med'])
    filtered = conns[mask]
    # remove recordings with gain errors
    mask = filtered['abs_deconv_base_amp_med'] < 0.02
    # remove recordings likely to have high crosstalk
    mask &= filtered['electrode_distance'] > 1
    # remove recordings with low sample count
    mask &= filtered['n_samples'] > 50


    filtered = filtered[mask]

    brushes = [pg.mkBrush('y') if c['synapse'] else pg.mkBrush(0.5) for c in filtered]

    c_mask = filtered['synapse']
    u_mask = ~c_mask

    signal = filtered['abs_deconv_amp_med']
    background = filtered['abs_deconv_base_amp_med']

    c_plt = scatter_plot.plot(background[c_mask], signal[c_mask], pen=None, symbol='d', symbolPen=None, symbolBrush='y', symbolSize=8)
    u_plt = scatter_plot.plot(background[u_mask], signal[u_mask], pen=None, symbol='o', symbolPen=None, symbolBrush=(200, 200, 200, 100), symbolSize=4)
    u_plt.setZValue(-10)
    # u_plt.scatter.setCompositionMode(pg.QtGui.QPainter.CompositionMode_Plus)

