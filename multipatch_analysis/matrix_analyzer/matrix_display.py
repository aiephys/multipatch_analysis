"""
Controls Matrix display and 2 side plots for main tab of Matrix Analyzer. 
Includes ColorMapping control for Matrix
Scatter plot displays histogram of matrix data as well as element data when element is selected
Trace plot displays average trace responses from each pair in a selected element

"""

from __future__ import print_function, division


import numpy as np
import pyqtgraph as pg
import pandas as pd
from analyzers import results_scatter, FormattableNumber
from multipatch_analysis.ui.graphics import MatrixItem
from pyqtgraph import parametertree as ptree
from pyqtgraph.parametertree import Parameter
from pyqtgraph.widgets.ColorMapWidget import ColorMapParameter

class SignalHandler(pg.QtCore.QObject):
        """Because we can't subclass from both QObject and QGraphicsRectItem at the same time
        """
        sigOutputChanged = pg.QtCore.Signal(object) #self

class MatrixTab(pg.QtGui.QWidget):
    def __init__(self):
        pg.QtGui.QWidget.__init__(self)
        self.layout = pg.QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.h_splitter = pg.QtGui.QSplitter()
        self.h_splitter.setOrientation(pg.QtCore.Qt.Horizontal)
        self.layout.addWidget(self.h_splitter, 0, 0)
        self.control_panel_splitter = pg.QtGui.QSplitter()
        self.control_panel_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.h_splitter.addWidget(self.control_panel_splitter)
        self.update_button = pg.QtGui.QPushButton("Update Results")
        self.control_panel_splitter.addWidget(self.update_button)
        self.ptree = ptree.ParameterTree(showHeader=False)
        self.control_panel_splitter.addWidget(self.ptree)
        self.matrix_widget = MatrixWidget()
        self.h_splitter.addWidget(self.matrix_widget)
        self.plot_splitter = pg.QtGui.QSplitter()
        self.plot_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.h_splitter.addWidget(self.plot_splitter)
        self.scatter_plot = MatrixScatterPlot()
        self.trace_plot = MatrixTracePlot()
        self.plot_splitter.addWidget(self.scatter_plot)
        self.plot_splitter.addWidget(self.trace_plot)
        self.h_splitter.setSizes([150, 300, 200])

class MatrixDisplayFilter(object):
    def __init__(self, view_box, data_fields, confidence_fields):
        self.output = None
        self._signalHandler = SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged
        self.view_box = view_box
        self.legend = None
        self.colorMap = ColorMapParameter()
        self.data_fields = data_fields
        self.confidence_fields = confidence_fields
        self.field_names = [field[0] for field in self.data_fields]

        self.colorMap.setFields(self.data_fields)

        self.params = Parameter.create(name='Matrix Display', type='group', children=[
            self.colorMap,
            {'name': 'Text format', 'type': 'str'},
            {'name': 'Show Confidence', 'type': 'list', 'values': self.confidence_fields},
            {'name': 'log_scale', 'type': 'bool'},
        ])
    
        self.params.sigTreeStateChanged.connect(self.invalidate_output)

    def colormap_legend(self):
        if self.legend is not None:
            self.view_box.removeItem(self.legend)
        cmap_item = [cmap for cmap in self.colorMap.children() if cmap['Enabled'] is True][0]
        log_scale = self.params.child('log_scale').value()
        colors = cmap_item.value().color
        x_min = cmap_item['Min']
        x_max = cmap_item['Max']
        x = np.linspace(x_min, x_max, len(colors))
        name = cmap_item.name()
        # units = self.colorMap.fields[name].get('units', None)
        scale, prefix = pg.siScale(x_min)
        # if units is not None:
        #     units = scale + units
        # else:
        #     units = ''
        self.legend = pg.GradientLegend([25, 300], [-20, -30])
        if log_scale is True:
            cmap2 = pg.ColorMap(x, colors)
            self.legend.setGradient(cmap2.getGradient())
            self.legend.setLabels({'%0.02f' % (a*scale):b for a,b in zip(cmap_item.value().pos, x)})
        else:
            self.legend.setGradient(cmap_item.value().getGradient())
            self.legend.setLabels({'%0.02f' % (a*scale):b for a,b in zip(x, cmap_item.value().pos)})
        self.view_box.addItem(self.legend)

    def element_display_output(self, result, mode='mean'):
        colormap = self.colorMap
        show_confidence = self.params['Show Confidence']
        text_format = self.params['Text format']
        

        if result[show_confidence] is not None:
            self.output = {'bordercolor': 0.6}
            default_bgcolor = np.array([128., 128., 128., 255.])
        else:
            self.output = {'bordercolor': 0.8}
            default_bgcolor = np.array([220., 220., 220.])
        
        if result['no_data'] is True:
            self.output['bgcolor'] = tuple(default_bgcolor)
            self.output['fgcolor'] = 0.6
            self.output['text'] = ''
        else:
            # mappable_result = {k:v for k,v in result.iteritems() if np.isscalar(v)}
            mappable_result = {}
            for k, v in result.iteritems():
                # d_type = np.array(v).dtype
                if np.isscalar(v):
                    mappable_result[k] = v
                elif isinstance(v, list) and k in self.field_names:
                    if mode=='mean':
                        mappable_result[k] = FormattableNumber(np.nanmean(v))
                else:
                    mappable_result[k] = np.nan

            color = colormap.map(mappable_result)[0]
        
            # desaturate low confidence cells
            if result[show_confidence] is not None:
                lower, upper = result[show_confidence]
                confidence = (1.0 - (upper - lower)) ** 2
                color = color * confidence + default_bgcolor * (1.0 - confidence)
            # invert text color for dark background
            self.output['fgcolor'] = 'w' if sum(color[:3]) < 384 else 'k'
            r2 = result.copy()
            r2.update(pd.Series(mappable_result))
            self.output['text'] = text_format.format(**r2)
            self.output['bgcolor'] = tuple(color)

        return self.output


    def invalidate_output(self):
        self.output = None


class MatrixWidget(pg.GraphicsLayoutWidget):
    sigClicked = pg.QtCore.Signal(object, object, object, object) # self, matrix_item, row, col
    def __init__(self):
        pg.GraphicsLayoutWidget.__init__(self)
        self.setRenderHints(self.renderHints() | pg.QtGui.QPainter.Antialiasing)
        v = self.addViewBox()
        v.setBackgroundColor('w')
        v.setAspectLocked()
        v.invertY()
        self.view_box = v
        self.matrix = None

    def set_matrix_data(self, text, fgcolor, bgcolor, border_color, rows, cols, size=50, header_color='k'):
        if self.matrix is not None:
            self.view_box.removeItem(self.matrix)

        self.matrix = MatrixItem(text=text, fgcolor=fgcolor, bgcolor=bgcolor, border_color=border_color,
                    rows=rows, cols=rows, size=50, header_color='k')
        self.matrix.sigClicked.connect(self.matrix_element_clicked)
        self.view_box.addItem(self.matrix)

    def matrix_element_clicked(self, matrix_item, event, row, col):
        self.sigClicked.emit(self, event, row, col) 


class MatrixScatterPlot(pg.GraphicsLayoutWidget):
    def __init__(self):
        pg.GraphicsLayoutWidget.__init__(self)
        self.setRenderHints(self.renderHints() | pg.QtGui.QPainter.Antialiasing)


class MatrixTracePlot(pg.GraphicsLayoutWidget):
    def __init__(self):
        pg.GraphicsLayoutWidget.__init__(self)
        self.setRenderHints(self.renderHints() | pg.QtGui.QPainter.Antialiasing)

class MatrixDisplay(object):
    sigClicked = pg.QtCore.Signal(object, object, object, object, object) # self, matrix_item, event, row, col
    def __init__(self, window, output_fields, confidence_fields, field_map, mode='mean'):
        self.matrix_tab = window
        self.matrix_widget = self.matrix_tab.matrix_widget
        self.matrix_display_filter = MatrixDisplayFilter(self.matrix_widget.view_box, output_fields, confidence_fields)
        self.field_map = field_map
        self.mode = mode
        self.matrix_map = {}
        self.hist_plot = None
        self.line = None
        self.scatter = None
        self.trace_plot = None
        self.trace_plot_list = []
        self.element = None
        self.selected = 0
        self.colors = [(0, 255, 0), (255, 0, 0), (0, 0, 255), (254, 169, 0), (170, 0, 127), (0, 230, 230)]

        self.matrix_widget.sigClicked.connect(self.display_matrix_element_data)

    def display_matrix_element_data(self, matrix_item, event, row, col):
        pre_class, post_class = self.matrix_map[row, col]
        analyzer = self.field_map[self.field_name]
        analyzer.print_element_info(pre_class, post_class, self.field_name)
        if self.hist_plot is not None:
            if int(event.modifiers() & pg.QtCore.Qt.ControlModifier)>0:
                self.selected += 1
                if self.selected >= len(self.colors):
                    self.selected = 0
                self.display_element_output(row, col, analyzer, trace_plot_list=self.trace_plot_list)
            else:
                self.display_element_reset() 
                self.display_element_output(row, col, analyzer)

    def display_element_output(self, row, col, analyzer, trace_plot_list=None):
        color = self.colors[self.selected]
        self.element = self.matrix_widget.matrix.cells[row][col]
        self.element.setPen(pg.mkPen({'color': color, 'width': 5}))
        pre_class, post_class = self.matrix_map[row, col]
        self.trace_plot = self.matrix_tab.trace_plot.addPlot()
        self.trace_plot_list.append(self.trace_plot)
        self.line, self.scatter = analyzer.plot_element_data(pre_class, post_class, self.field_name, color=color, trace_plt=self.trace_plot)
        if len(self.trace_plot_list) > 1:
            first_plot = self.trace_plot_list[0]
            for plot in self.trace_plot_list[1:]:
                plot.setYLink(first_plot)
                plot.setXLink(first_plot)
        self.hist_plot.addItem(self.line)
        if self.scatter is not None:
            self.hist_plot.addItem(self.scatter)

    def display_element_reset(self):
        self.selected = 0
        if self.hist_plot is not None:
            [self.hist_plot.removeItem(item) for item in self.hist_plot.items[1:]]
        if self.trace_plot is not None:
           self.matrix_tab.trace_plot.clear()
           self.trace_plot = None
        self.trace_plot_list = []

        show_confidence = self.matrix_display_filter.params['Show Confidence']
        bordercolor = 0.6 if show_confidence is None else 0.8
        for cells in self.matrix_tab.matrix_widget.matrix.cells:
            for cell in cells:
                cell.setPen(pg.mkPen({'color': bordercolor, 'width': 1}))

    def update_matrix_display(self, results, cell_classes, cell_groups):
        
        shape = (len(cell_groups),) * 2
        text = np.empty(shape, dtype=object)
        fgcolor = np.empty(shape, dtype=object)
        bgcolor = np.empty(shape, dtype=object)
        bordercolor = np.empty(shape, dtype=object)
        self.matrix_display_filter.colormap_legend()

        # call display function on every matrix element
        
        for i,row in enumerate(cell_groups):
            for j,col in enumerate(cell_groups):
                # output = self.matrix_display_filter.element_display_output(results[(row, col)])
                output = self.matrix_display_filter.element_display_output(results[(row, col)], mode=self.mode)
                self.matrix_map[i, j] = (row, col)
                text[i, j] = output['text']
                fgcolor[i, j] = output['fgcolor']
                bgcolor[i, j] = output['bgcolor']
                bordercolor[i, j] = output['bordercolor']
                
        # Force cell class descriptions down to tuples of 2 items
        # Kludgy, but works for now.
        # update 3/8/19: Doesn't work for CellClasses of 1 item,
        # attempt to fix so it doesn't break in mp_a\ui\graphics.py
        # at line 90. 
        rows = []
        for i,cell_class in enumerate(cell_classes):
            tup = cell_class.as_tuple
            row = tup[:1]
            if len(tup) > 1:
                row = row + (' '.join(tup[1:]),)
            else:
                row = (' '*i,) + row
            # if len(tup) > 1:
            #     row = tup
            # elif len(tup) == 1:
            #     row = list(tup)
            rows.append(row)

        self.matrix_tab.matrix_widget.set_matrix_data(text=text, fgcolor=fgcolor, bgcolor=bgcolor, border_color=bordercolor,
                    rows=rows, cols=rows, size=50, header_color='k')


        # plot hist or scatter of data in side panel, if there are multiple colormaps take the first enabled one
        color_map_fields = self.matrix_display_filter.colorMap.children()
        if len(color_map_fields) > 1:
            self.field_name = [field.name() for field in color_map_fields if field['Enabled'] is True][0]
        else:
            self.field_name = color_map_fields[0].name()
        field = self.matrix_display_filter.colorMap.fields[self.field_name]
        if self.hist_plot is not None:
            self.matrix_tab.scatter_plot.removeItem(self.hist_plot)
        self.hist_plot = self.matrix_tab.scatter_plot.addPlot()
        vals = results.loc[self.field_name, :]
        vals = seriesToArray(vals)
        y, x = np.histogram(vals, bins=np.linspace(min(vals), max(vals), 10))
        self.hist_plot.plot(x, y, stepMode=True, fillLevel=0, brush=(255,255,255,150))
        units = field.get('units', '')
        self.hist_plot.setLabels(left='Count', bottom=(self.field_name, units))

        # self.analysis.summary(self.results, self.field_name)

def seriesToArray(series):
    ## convert pandas series to array
    s = series[series.notnull()]
    s_list = s.to_list()
    s_array = [y for x in s_list if isinstance(x, list) for y in x]
    if len(s_array) == 0:
        return s
    s_array = [val for val in s_array if ~np.isnan(val)]
    return s_array