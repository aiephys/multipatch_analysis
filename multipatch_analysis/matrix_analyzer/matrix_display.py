"""
Controls Matrix display. 
Includes ColorMapping control for Matrix

"""

from __future__ import print_function, division


import numpy as np
import pyqtgraph as pg
import pandas as pd
from collections import OrderedDict
from analyzers import results_scatter, FormattableNumber
from multipatch_analysis.ui.graphics import MatrixItem
from pyqtgraph import parametertree as ptree
from pyqtgraph.parametertree import Parameter
from pyqtgraph.widgets.ColorMapWidget import ColorMapParameter

class SignalHandler(pg.QtCore.QObject):
        """Because we can't subclass from both QObject and QGraphicsRectItem at the same time
        """
        sigOutputChanged = pg.QtCore.Signal(object) #self



class MatrixDisplayFilter(object):
    def __init__(self, view_box, data_fields):
        self.output = None
        self._signalHandler = SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged
        self.view_box = view_box
        self.legend = None
        self.colorMap = ColorMapParameter()
        self.data_fields = data_fields
        self.field_names = [field[0] for field in self.data_fields]

        self.colorMap.setFields(self.data_fields)

        self.params = Parameter.create(name='Matrix Display', type='group', children=[
            self.colorMap,
            {'name': 'Text format', 'type': 'str'},
            {'name': 'Show Confidence', 'type': 'list', 'values': [field[0] for field in self.data_fields]},
            {'name': 'log_scale', 'type': 'bool'},
        ])
    
        self.params.sigTreeStateChanged.connect(self.invalidate_output)

    def element_display_output(self, result, default_bgcolor):
        colormap = self.colorMap
        show_confidence = self.params['Show Confidence']
        text_format = self.params['Text format']
        self.output = {}

        # if show_confidence != 'None':
        #     default_bgcolor = np.array([128., 128., 128., 255.])
        # else:
        #     default_bgcolor = np.array([220., 220., 220.])
        
        # if result['no_data'] is False:
        #     self.output['bgcolor'] = tuple(default_bgcolor)
        #     self.output['fgcolor'] = 0.6
        #     self.output['text'] = ''
        # else:
        result_vals = result.loc[:, 'metric_summary']
        mappable_result = {k:v for k,v in result_vals.iteritems() if np.isscalar(v)}

        color = colormap.map(mappable_result)[0]
    
        # desaturate low confidence cells
        if show_confidence != 'None':
            lower, upper = result[show_confidence, 'metric_conf']
            confidence = (1.0 - (upper - lower)) ** 2
            color = color * confidence + default_bgcolor * (1.0 - confidence)
        # invert text color for dark background
        self.output['fgcolor'] = 'w' if sum(color[:3]) < 384 else 'k'
        text_result = {k:FormattableNumber(v) if isinstance(v, float) else v for k, v in result_vals.iteritems()}
        self.output['text'] = text_format.format(**text_result)
        self.output['bgcolor'] = tuple(color)

        return self.output
    
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

    def get_colormap_field(self):
        color_map_fields = self.colorMap.children()
        if len(color_map_fields) > 1:
            field_name = [field.name() for field in color_map_fields if field['Enabled'] is True][0]
        else:
            field_name = color_map_fields[0].name()

        return field_name

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
                    rows=rows, cols=cols, size=50, header_color='k')
        self.matrix.sigClicked.connect(self.matrix_element_clicked)
        self.view_box.addItem(self.matrix)

    def matrix_element_clicked(self, matrix_item, event, row, col):
        self.sigClicked.emit(self, event, row, col) 


class MatrixDisplay(object):
    sigClicked = pg.QtCore.Signal(object, object, object, object, object) # self, matrix_item, event, row, col
    def __init__(self, window, output_fields, field_map):
        self.main_window = window
        self.matrix_widget = self.main_window.matrix_widget
        self.matrix_display_filter = MatrixDisplayFilter(self.matrix_widget.view_box, output_fields)
        self.field_map = field_map
        self.element = None

    def get_element_classes(self, row, col):
        pre_class, post_class = [k for k, v in self.matrix_map.items() if v==[row, col]][0]
        return pre_class, post_class

    def get_field_analyzer(self):
        analyzer = self.field_map[self.field_name]
        return analyzer, field_name

    def color_element(self, row, col, color):
        self.element = self.matrix_widget.matrix.cells[row[0]][col[0]]
        self.element.setPen(pg.mkPen({'color': color, 'width': 5}))

    def element_color_reset(self):
        show_confidence = self.matrix_display_filter.params['Show Confidence']
        bordercolor = 0.6 if show_confidence is None else 0.8
        for cells in self.main_window.matrix_widget.matrix.cells:
            for cell in cells:
                cell.setPen(pg.mkPen({'color': bordercolor, 'width': 1}))

    def update_matrix_display(self, results, group_results, cell_classes, cell_groups, field_map, pre_cell_classes=None, post_cell_classes=None):
        self.results = results
        self.group_results = group_results
        self.matrix_map = OrderedDict()
        show_confidence = self.matrix_display_filter.params['Show Confidence']

        shape = (len(cell_groups),) * 2
        text = np.empty(shape, dtype=object)
        text.fill(u'')
        fgcolor = np.empty(shape, dtype=object)
        fgcolor.fill(0.6)
        bgcolor = np.empty(shape, dtype=object)
        if show_confidence != 'None':
            default_bordercolor = 0.6
            # default_bgcolor = np.array([128., 128., 128., 255.])
        else:
            default_bordercolor = 0.8
            # default_bgcolor = np.array([220., 220., 220.])
        default_bgcolor = np.array([128., 128., 128., 255.])
        bgcolor.fill(tuple(default_bgcolor))
        bordercolor = np.empty(shape, dtype=object)
        bordercolor.fill(default_bordercolor)
        self.matrix_display_filter.colormap_legend()

        ## for backwards compatability and to avoid an api change, make it okay to just pass cell_classes
        if pre_cell_classes == None:
            pre_cell_classes = cell_classes
        if post_cell_classes == None:
            post_cell_classes = cell_classes

        # call display function on every matrix element
        
        for ii, pre in enumerate(pre_cell_classes):
            for jj, post in enumerate(post_cell_classes):
                self.matrix_map[(pre, post)] = [ii, jj]
        for group, result in self.group_results.iterrows():
            try:
                ## do some shenanigans because cell_class == cell_class doens't work, instead need cell_class.name == cell_class.name
                group_name = (group[0].name, group[1].name)
                keys_names = [(pre.name, post.name) for pre, post in self.matrix_map.keys()]
                key_index = keys_names.index(group_name)
                i, j = self.matrix_map[self.matrix_map.keys()[key_index]]
            except ValueError:
                ### some group_results catagories won't be in matrix_map because of pre and post cell class specifications
                continue

            no_data = all([result.get('conn_no_data',{}).get('metric_summary', True), result.get('strength_no_data',{}).get('metric_summary', True), result.get('dynamics_no_data',{}).get('metric_summary', True)])
            if no_data is False:
                output = self.matrix_display_filter.element_display_output(result, default_bgcolor)
                text[i, j] = output['text']
                fgcolor[i, j] = output['fgcolor']
                bgcolor[i, j] = output['bgcolor']
                
        # Force cell class descriptions down to tuples of 2 items
        # Kludgy, but works for now.
        # update 3/8/19: Doesn't work for CellClasses of 1 item,
        # attempt to fix so it doesn't break in mp_a\ui\graphics.py
        # at line 90. 
        rows = []
        cols = []
        pre_names = [pre.name for pre in pre_cell_classes] ## equality comparison happens on names, so need to convert to list of names here to ask if cell is in the list later
        post_names = [post.name for post in post_cell_classes]
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
            if cell_class in pre_names:
                rows.append(row)
            if cell_class in post_names:
                cols.append(row)

        self.main_window.matrix_widget.set_matrix_data(text=text, fgcolor=fgcolor, bgcolor=bgcolor, border_color=bordercolor,
                    rows=rows, cols=cols, size=50, header_color='k')

        # self.analysis.summary(self.results, self.field_name)
