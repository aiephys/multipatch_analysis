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

    def set_fields(self,fields, field_map):
        self.fields = fields
        self.field_map = field_map
        self.setFields(self.fields)

    def set_data(self, data):
        field_names = self.fields.keys()
        # fm = {k:v.name for k, v in self.field_map.items()}
        # fm_df = pd.DataFrame.from_dict(fm, orient='index').loc[field_names]
        # fm_df.columns = ['analyzer']
        # field_data = data.join(fm_df, how='inner')
        # field_data.drop(columns=['analyzer'], inplace=True)
        field_data = data.reindex(field_names)
        summary_data = field_data.applymap(lambda x: result_mean(x))
        rec_data = summary_data.transpose().to_records()
        names = ('pre_class', 'post_class') + tuple([str(name) for name in rec_data.dtype.names[2:]])
        rec_data.dtype.names = names

        self.setData(rec_data)

    def invalidate_output(self):
        self.data = None

class PairScatterPlot(pg.ScatterPlotWidget):
    def __init__(self):
        pg.ScatterPlotWidget.__init__(self)

        header = pg.QtGui.QLabel()
        header.setText('<span style="font-weight: bold">Pair-wise Scatter Plot</span>')
        self.ctrlPanel.insertWidget(0, header)

    def set_fields(self, fields, field_map):
        self.fields = fields
        self.field_map = field_map
        self.setFields(fields)

    def set_data(self, data):
        field_names = self.fields.keys()
        field_data = data.reindex(field_names)
        rec_data = field_data.transpose().to_records()
        names = ('pre_class', 'post_class') + tuple([str(name) for name in rec_data.dtype.names[2:]])
        rec_data.dtype.names = names

        self.setData(rec_data)
        
    def invalidate_output(self):
        self.data = None

def result_mean(result):
    if np.isscalar(result):
        return result
    else:
        return np.nanmean(result)

