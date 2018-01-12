# *-* coding: utf-8 *-*
from __future__ import print_function, division
import numpy as np
import pyqtgraph as pg
from neuroanalysis.stats import binomial_ci
from neuroanalysis.ui.plot_grid import PlotGrid


class MatrixItem(pg.QtGui.QGraphicsItemGroup):
    """GraphicsItem displaying a table with column / row labels and text in
    each cell.
    
    Parameters
    ----------
    text : 2d array or nested lists
        Strings to display inside each cell
    fgcolor : 2d array or nested lists
        Text colors for each cell
    bgcolor : 2d array or nested lists
        Background colors for each cell
    rows : 1d array or list
        Strings to display as row header
    cols : 1d array or list
        Strings to display as col header
    size : float
        Width of each cell
    header_color : str | tuple
        Color of header text
    """
    def __init__(self, text, fgcolor, bgcolor, rows, cols, size=50, header_color='w'):
        pg.QtGui.QGraphicsItemGroup.__init__(self)

        for i,row in enumerate(rows):
            txt = pg.QtGui.QGraphicsTextItem(str(row), parent=self)
            br = txt.boundingRect()
            txt.setPos(-br.width() - 10, i * size + size/2. - br.center().y())
            txt.setDefaultTextColor(pg.mkColor(header_color))

        for i,col in enumerate(cols):
            txt = pg.QtGui.QGraphicsTextItem(str(col), parent=self)
            txt.rotate(90)
            br = txt.mapRectToParent(txt.boundingRect())
            txt.setPos(i * size + size/2 - br.center().x(), -br.height() - 10)
            txt.setDefaultTextColor(pg.mkColor(header_color))

        for i,row in enumerate(rows):
            for j,col in enumerate(cols):
                x = j * size
                y = i * size
                rect = pg.QtGui.QGraphicsRectItem(x, y, size, size, parent=self)
                rect.setBrush(pg.mkBrush(bgcolor[i][j]))
                rect.setZValue(-10)

                txt = pg.QtGui.QGraphicsTextItem(text[i][j], parent=self)
                br = txt.boundingRect()
                txt.setPos(x + size/2 - br.center().x(), y + size/2 - br.center().y())
                txt.setDefaultTextColor(pg.mkColor(fgcolor[i][j]))

        br = pg.QtCore.QRectF()
        for item in self.childItems():
            br = br.united(self.mapRectFromItem(item, item.boundingRect()))
        self._bounding_rect = br

    def boundingRect(self):
        return self._bounding_rect
    
    
def distance_plot(connected, distance, plots=None, color=(100, 100, 255), window=40e-6, spacing=None, name=None, fill_alpha=30):
    """Draw connectivity vs distance profiles with confidence intervals.
    
    Parameters
    ----------
    connected : boolean array
        Whether a synaptic connection was found for each probe
    distance : array
        Distance between cells for each probe
    plots : list of PlotWidget | PlotItem
        (optional) Two plots used to display distance profile and scatter plot.
    """
    color = pg.colorTuple(pg.mkColor(color))[:3]
    connected = np.array(connected).astype(float)
    distance = np.array(distance)
    pts = np.vstack([distance, connected]).T
    
    # scatter points a bit
    conn = pts[:,1] == 1
    unconn = pts[:,1] == 0
    if np.any(conn):
        cscat = pg.pseudoScatter(pts[:,0][conn], spacing=10e-6, bidir=False)
        mx = abs(cscat).max()
        if mx != 0:
            cscat = cscat * 0.2# / mx
        pts[:,1][conn] = -2e-5 - cscat
    if np.any(unconn):
        uscat = pg.pseudoScatter(pts[:,0][unconn], spacing=10e-6, bidir=False)
        mx = abs(uscat).max()
        if mx != 0:
            uscat = uscat * 0.2# / mx
        pts[:,1][unconn] = uscat

    # scatter plot connections probed
    if plots is None:
        grid = PlotGrid()
        grid.set_shape(2, 1)
        grid.grid.ci.layout.setRowStretchFactor(0, 5)
        grid.grid.ci.layout.setRowStretchFactor(1, 10)
        plots = (grid[1,0], grid[0,0])
        plots[0].grid = grid
        plots[0].addLegend()
        grid.show()
    plots[0].setLabels(bottom=('distance', 'm'), left='connection probability')

    if plots[1] is not None:
        plots[1].setXLink(plots[0])
        plots[1].hideAxis('bottom')
        plots[1].hideAxis('left')

        color2 = color + (100,)
        scatter = plots[1].plot(pts[:,0], pts[:,1], pen=None, symbol='o', labels={'bottom': ('distance', 'm')}, symbolBrush=color2, symbolPen=None, name=name)
        scatter.scatter.opts['compositionMode'] = pg.QtGui.QPainter.CompositionMode_Plus

    # use a sliding window to plot the proportion of connections found along with a 95% confidence interval
    # for connection probability

    if spacing is None:
        spacing = window / 4.0
        
    xvals = np.arange(window / 2.0, 500e-6, spacing)
    upper = []
    lower = []
    prop = []
    ci_xvals = []
    for x in xvals:
        minx = x - window / 2.0
        maxx = x + window / 2.0
        # select points inside this window
        mask = (distance >= minx) & (distance <= maxx)
        pts_in_window = connected[mask]
        # compute stats for window
        n_probed = pts_in_window.shape[0]
        n_conn = pts_in_window.sum()
        if n_probed == 0:
            prop.append(np.nan)
        else:
            prop.append(n_conn / n_probed)
            ci = binomial_ci(n_conn, n_probed)
            lower.append(ci[0])
            upper.append(ci[1])
            ci_xvals.append(x)

    # plot connection probability and confidence intervals
    color2 = [c / 3.0 for c in color]
    mid_curve = plots[0].plot(xvals, prop, pen={'color': color, 'width': 3}, antialias=True, name=name)
    upper_curve = plots[0].plot(ci_xvals, upper, pen=(0, 0, 0, 0), antialias=True)
    lower_curve = plots[0].plot(ci_xvals, lower, pen=(0, 0, 0, 0), antialias=True)
    upper_curve.setVisible(False)
    lower_curve.setVisible(False)
    color2 = color + (fill_alpha,)
    fill = pg.FillBetweenItem(upper_curve, lower_curve, brush=color2)
    fill.setZValue(-10)
    plots[0].addItem(fill, ignoreBounds=True)
    
    return plots