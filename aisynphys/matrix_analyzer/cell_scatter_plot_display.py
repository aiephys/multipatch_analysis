"""
ScatterPlotWidgets for Matrix Analyzer. 
Allows plotting of Matrix data on a per cell basis.

"""

from __future__ import print_function, division
import pyqtgraph as pg
import colorsys
import pandas as pd
import numpy as np

class CellScatterTab(pg.QtGui.QWidget):
    def __init__(self):
        pg.QtGui.QWidget.__init__(self)
        self.layout = pg.QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.v_splitter = pg.QtGui.QSplitter()
        # self.v_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        # self.layout.addWidget(self.v_splitter)
        self.cell_scatter = CellScatterPlot()
        self.layout.addWidget(self.cell_scatter)
        # self.morpho_scatter = MorphoScatter()
        # self.v_splitter.add(self.morpho_scatter)

        # self.colorMap.sigColorMapChanged.connect(self.set_colors)

class CellScatterPlot(pg.ScatterPlotWidget):
    def __init__(self):
        pg.ScatterPlotWidget.__init__(self)

    def set_fields(self, fields):
        self.fields = [('CellClass', {'mode': 'enum'})]
        self.fields.extend([f for f in fields if f != ('None', {})])
        self.setFields(self.fields)

    def set_data(self, data):
        rec_data = data.to_records()
        names = tuple([str(name) for name in rec_data.dtype.names]) # for some reason the unicode throws off
        rec_data.dtype.names = names

        for field, options in self.fields.items():
            data_type = options.get('mode')
            values = options.get('values')
            defaults = options.get('defaults')
            if data_type == 'enum':
                if values is None:
                    self.fields[field]['values'] = list(set(data.get(field)))
                if defaults is None:
                    n_colors = len(set(data.get(field))) if values is None else len(values)
                    self.fields[field]['defaults'] = {'colormap': [pg.intColor(n, n_colors) for n in np.arange(n_colors)]}
                
        self.setData(rec_data)

    # def map(self, data):
    #     if self.mapType == 'range':
    #         pg.ColorMap

    def invalidate_output(self):
        self.data = None

# class PatchSeqScatter(CellScatterPlot):
#     def __init__(self):
#         CellScatterPlot.__init__(self)

# class MorphoScatter(CellScatterPlot):
#     def __init__(self):
#         CellScatterPlot.__init__(self)