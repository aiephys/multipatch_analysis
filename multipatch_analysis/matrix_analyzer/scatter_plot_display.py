"""
ScatterPlotWidgets for Matrix Analyzer. 
Allows plotting of Matrix data on a per element basis and per pair basis.

"""

from __future__ import print_function, division
import pyqtgraph as pg

class ScatterPlotTab(pg.QtGui.QWidget):
    def __init__(self):
        pg.QtGui.QWidget.__init__(self)
        self.layout = pg.QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.h_splitter = pg.QtGui.QSplitter()
        self.h_splitter.setOrientation(pg.QtCore.Qt.Horizontal)
        self.layout.addWidget(self.h_splitter)
        self.element_scatter = ElementScatterPlot()
        self.pair_scatter = PairScatterPlot()
        self.h_splitter.addWidget(self.element_scatter)
        self.h_splitter.addWidget(self.pair_scatter)

class ElementScatterPlot(pg.ScatterPlotWidget):
    def __init__(self):
        pg.ScatterPlotWidget.__init__(self)

        header = pg.QtGui.QLabel()
        header.setText('<span style="font-weight: bold">Element-wise Scatter Plot</span>')
        self.ctrlPanel.insertWidget(0, header)

    def set_fields(self, fields):
        self.setFields(fields)

    def set_data(self, data):
        self.setData(data)

    def invalidate_output(self):
        self.data = None

class PairScatterPlot(pg.ScatterPlotWidget):
    def __init__(self):
        pg.ScatterPlotWidget.__init__(self)

        header = pg.QtGui.QLabel()
        header.setText('<span style="font-weight: bold">Pair-wise Scatter Plot</span>')
        self.ctrlPanel.insertWidget(0, header)

    def set_fields(self, fields):
        self.setFields(fields)

    def set_data(self, data):
        self.setData(data)
        
    def invalidate_output(self):
        self.data = None