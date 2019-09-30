from __future__ import division
import numpy as np
import pyqtgraph as pg
from aisynphys.experiment_list import cached_experiments

expts = cached_experiments()


dates = [(e.date - expts[0].date).days for e in expts]
ud = list(np.unique(dates))
days = [ud.index(d) for d in dates]

probed = np.array([e.n_connections_probed for e in expts])
found = np.array([e.n_connections for e in expts])

ctypes = [tuple(exp.cre_types) for exp in expts]
n_ctypes = len(set(ctypes))
ctype_colors = {}
for exp in expts:
    ct = tuple(exp.cre_types)
    color = pg.intColor(len(ctype_colors), hues=6, values=2, minValue=100)
    ctype_colors.setdefault(ct, pg.mkBrush(color))

expt_colors = [ctype_colors[ct] for ct in ctypes]


p1 = pg.plot(dates, np.cumsum(probed), symbol='o', symbolBrush=expt_colors, symbolPen=None, title="Connections probed over time (including no-experiment days)")

p1 = pg.plot(days, np.cumsum(probed), symbol='o', symbolBrush=expt_colors, symbolPen=None, title="Connections probed per experiment day")

p2 = pg.plot(dates, np.cumsum(found), symbol='o', symbolBrush=expt_colors, symbolPen=None, title="Synapses detected over time")
p2.plot(days, np.cumsum(found), symbol='o', symbolBrush=expt_colors, symbolPen=None)

p3 = pg.plot(dates, found/probed, pen=None, symbol='o', symbolBrush=expt_colors, symbolPen=None, title="Connectivity over time")


from pyqtgraph.graphicsItems.LegendItem import ItemSample
import pyqtgraph.functions as fn
from pyqtgraph.graphicsItems.ScatterPlotItem import drawSymbol
from pyqtgraph.Qt import QtGui, QtCore


class Item(ItemSample):
    def __init__(self, **kwds):
        self.opts = kwds
        ItemSample.__init__(self, self)

    def boundingRect(self):
        return QtCore.QRectF(0, 0, 20, 20)
        
    def paint(self, p, *args):
        opts = self.opts
        if opts.get('antialias', True):
            p.setRenderHint(p.Antialiasing)
        
        if opts.get('fillLevel', None) is not None and opts.get('fillBrush', None) is not None:
            p.setBrush(fn.mkBrush(opts['fillBrush']))
            p.setPen(fn.mkPen(None))
            p.drawPolygon(QtGui.QPolygonF([
                QtCore.QPointF(2, 18), 
                QtCore.QPointF(18, 2), 
                QtCore.QPointF(18, 18)]))
        
        if opts.get('pen', None) is not None:
            p.setPen(fn.mkPen(opts['pen']))
            p.drawLine(2, 18, 18, 2)
        
        symbol = opts.get('symbol', None)
        if symbol is not None:
                
            pen = fn.mkPen(opts.get('symbolPen', None))
            brush = fn.mkBrush(opts.get('symbolBrush', None))
            size = opts.get('symbolSize', 10)
            
            p.translate(10,10)
            path = drawSymbol(p, symbol, size, pen, brush)




legend = p1.addLegend()
for ct, color in ctype_colors.items():
    legend.addItem(Item(pen='w', symbol='o', symbolPen=None, symbolBrush=color, brush=None, size=10), '/'.join(ct))
