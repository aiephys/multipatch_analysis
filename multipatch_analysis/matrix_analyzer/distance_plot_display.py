"""
Distance vs connection probability plots for Matrix Analyzer.

"""

from __future__ import print_function, division
import pyqtgraph as pg
import numpy as np
from pyqtgraph.widgets.ColorMapWidget import ColorMapParameter
from pyqtgraph import parametertree as ptree
from neuroanalysis.ui.plot_grid import PlotGrid


class DistancePlotTab(pg.QtGui.QWidget):
    def __init__(self):
        pg.QtGui.QWidget.__init__(self)
        self.layout = pg.QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.h_splitter = pg.QtGui.QSplitter()
        self.h_splitter.setOrientation(pg.QtCore.Qt.Horizontal)
        self.layout.addWidget(self.h_splitter)
        self.distance_filter = DistanceFilter()
        self.ptree = ptree.ParameterTree(showHeader=False)
        self.ptree.setParameters(self.distance_filter.cmap)
        self.h_splitter.addWidget(self.ptree)
        self.distance_plot = DistancePlot()
        self.h_splitter.addWidget(self.distance_plot.grid)


class DistancePlot(object):
    def __init__(self):
        self.grid = PlotGrid()
        self.grid.set_shape(2, 1)
        self.grid.grid.ci.layout.setRowStretchFactor(0, 5)
        self.grid.grid.ci.layout.setRowStretchFactor(1, 10)
        self.plots = (self.grid[1,0], self.grid[0,0])
        self.plots[0].grid = self.grid
        self.plots[0].addLegend()
        self.grid.show()
        self.plots[0].setLabels(bottom=('distance', 'm'), left='connection probability')     

class DistanceFilter(object):
    def __init__(self):
        self.cmap = ColorMapParameter()
        
    def set_fields(self, pre_classes, post_classes):
        self.pre_class_list = pre_classes
        self.post_class_list = post_classes
        self.fields = [
            ('pre_class', {'mode': 'enum', 'values': self.pre_class_list}),
            ('post_class', {'mode': 'enum', 'values': self.post_class_list}),
        ]

        self.cmap.setFields(fields)

