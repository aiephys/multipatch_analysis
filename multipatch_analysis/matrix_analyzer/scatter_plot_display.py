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
        self.pair_scatter = PairScatterPlot()#self.element_scatter.plot)
        self.v_splitter.addWidget(self.element_scatter)
        self.v_splitter.addWidget(self.pair_scatter)

class ScatterPlots(pg.ScatterPlotWidget):
    def __init__(self):
        pg.ScatterPlotWidget.__init__(self)

    def set_fields(self, fields):
        self.fields = [('pair_class', {'mode': 'enum'})]
        self.fields.extend([f for f in fields if f != ('None', {})])
        self.setFields(self.fields)

    def color_selected_element(self, color, pre_class, post_class):
        pair_name = '-'.join([pre_class.name, post_class.name])
        try:
            pair_map = self.colorMap.child('pair_class')
        except KeyError:
            pair_map = self.colorMap.addNew('pair_class')
        pair_map['Values', pair_name] = pg.mkColor(color)

    def reset_element_color(self):
        try:
            pair_class_map = self.colorMap.child('pair_class')
            self.colorMap.removeChild(pair_class_map)
        except:
            return

    def invalidate_output(self):
        self.data = None

class ElementScatterPlot(ScatterPlots):
    def __init__(self):
        ScatterPlots.__init__(self)
    
        header = pg.QtGui.QLabel()
        header.setText('<span style="font-weight: bold">Element-wise Scatter Plot</span>')
        self.ctrlPanel.insertWidget(0, header)

    def set_data(self, data):
        field_data = data.xs('metric_summary', axis='columns', level=1, drop_level=True)
        field_data.reset_index(inplace=True)
        field_data['pair_class'] = field_data.apply(lambda row: '-'.join([row.pre_class.name, row.post_class.name]), axis=1)
        field_data.drop(columns=['pre_class', 'post_class'])
        rec_data = field_data.to_records()
        names = tuple([str(name) for name in rec_data.dtype.names]) # for some reason the unicode throws off
        rec_data.dtype.names = names

        self.fields['pair_class']['values'] = list(field_data.pair_class)
        self.setData(rec_data)

class PairScatterPlot(ScatterPlots):
    def __init__(self):
        ScatterPlots.__init__(self)
        
        header = pg.QtGui.QLabel()
        header.setText('<span style="font-weight: bold">Pair-wise Scatter Plot</span>')
        self.ctrlPanel.insertWidget(0, header)

    def set_data(self, data):
        data['pair_class'] = data.apply(lambda row: '-'.join([row.pre_class.name, row.post_class.name]), axis=1)
        rec_data = data.to_records()
        names = tuple([str(name) for name in rec_data.dtype.names]) # for some reason the unicode throws off
        rec_data.dtype.names = names

        self.fields['pair_class']['values'] = list(data.pair_class)
        self.setData(rec_data)