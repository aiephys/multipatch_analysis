"""
ScatterPlotWidgets for Matrix Analyzer. 
Allows plotting of Matrix data on a per element basis and per pair basis.

"""

from __future__ import print_function, division
import pyqtgraph as pg
import pandas as pd
import numpy as np

class ScatterPlotTab(pg.QtGui.QWidget):
    def __init__(self):
        pg.QtGui.QWidget.__init__(self)
        self.layout = pg.QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.v_splitter = pg.QtGui.QSplitter()
        self.v_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.layout.addWidget(self.v_splitter)
        self.element_scatter = ElementScatterPlot()
        self.pair_scatter = PairScatterPlot(plot=None)#self.element_scatter.plot)
        self.v_splitter.addWidget(self.element_scatter)
        self.v_splitter.addWidget(self.pair_scatter)

class ElementScatterPlot(pg.ScatterPlotWidget):
    def __init__(self):
        pg.ScatterPlotWidget.__init__(self)
        

        header = pg.QtGui.QLabel()
        header.setText('<span style="font-weight: bold">Element-wise Scatter Plot</span>')
        self.ctrlPanel.insertWidget(0, header)

    def set_fields(self, fields):
        self.fields = [f for f in fields if f != ('None', {})]
        self.setFields(self.fields)

    def set_data(self, data):
        field_data = data.xs('metric_summary', axis='columns', level=1, drop_level=True)
        field_data.reset_index(in_place=True)
        field_data['pair_class'] = list(zip(field_data.pre_class, field_data.post_class))
        field_data.drop(columns=['pre_class', 'post_class'])
        asdf
        rec_data = field_data.to_records()
        names = tuple([str(name) for name in rec_data.dtype.names]) # for some reason the unicode throws off
        rec_data.dtype.names = names

        self.setData(rec_data)

    def invalidate_output(self):
        self.data = None

class PairScatterPlot(pg.ScatterPlotWidget):
    def __init__(self, plot):
        pg.ScatterPlotWidget.__init__(self, plot=plot)

        header = pg.QtGui.QLabel()
        header.setText('<span style="font-weight: bold">Pair-wise Scatter Plot</span>')
        self.ctrlPanel.insertWidget(0, header)

    def set_fields(self, fields):
        self.fields = [f for f in fields if f != ('None', {})]
        self.setFields(self.fields)

    def set_data(self, data):
        rec_data = data.to_records()
        names = tuple([str(name) for name in rec_data.dtype.names]) # for some reason the unicode throws off
        rec_data.dtype.names = names

        self.setData(rec_data)
        
    def invalidate_output(self):
        self.data = None