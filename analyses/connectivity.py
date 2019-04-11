# -*- coding: utf-8 -*-

"""
Prototype code for analyzing connectivity and synaptic properties between cell classes.


"""

from __future__ import print_function, division

from collections import OrderedDict
import numpy as np
import pyqtgraph as pg
import multipatch_analysis.database as db
import re
from multipatch_analysis.connectivity import query_pairs, ConnectivityAnalyzer, StrengthAnalyzer, DynamicsAnalyzer, results_scatter, get_all_output_fields
from multipatch_analysis import constants
from multipatch_analysis.cell_class import CellClass, classify_cells, classify_pairs
from multipatch_analysis.ui.graphics import MatrixItem
from pyqtgraph import parametertree as ptree
from pyqtgraph.parametertree import Parameter
from pyqtgraph.widgets.ColorMapWidget import ColorMapParameter
from pyqtgraph.widgets.DataFilterWidget import DataFilterParameter


class MainWindow(pg.QtGui.QWidget):
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
        self.update_button = pg.QtGui.QPushButton("Update Matrix")
        self.control_panel_splitter.addWidget(self.update_button)
        self.ptree = ptree.ParameterTree(showHeader=False)
        self.control_panel_splitter.addWidget(self.ptree)
        # self.dsp = DataScatterPlot()
        # header = pg.QtGui.QLabel()
        # header.setText('<span style="font-weight: bold">Scatter Plot Display</span>')
        # self.control_panel_splitter.addWidget(header)
        # self.control_panel_splitter.addWidget(self.dsp.spw.ctrlPanel)
        # self.control_panel_splitter.setSizes([10, 300, 5, 200])

        self.matrix_widget = MatrixWidget()
        self.h_splitter.addWidget(self.matrix_widget)
        self.plot_splitter = pg.QtGui.QSplitter()
        self.plot_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.h_splitter.addWidget(self.plot_splitter)
        self.scatter_plot = ScatterPlot()
        self.trace_plot = TracePlot()
        # self.trace_plot.setVisible(False)
        self.plot_splitter.addWidget(self.scatter_plot)
        self.plot_splitter.addWidget(self.trace_plot)
        self.h_splitter.setSizes([150, 300, 200])


class SignalHandler(pg.QtCore.QObject):
        """Because we can't subclass from both QObject and QGraphicsRectItem at the same time
        """
        sigOutputChanged = pg.QtCore.Signal(object) #self

class ExperimentFilter(object):
    def __init__(self):  
        s = db.Session()
        self._signalHandler = SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged
        self.pairs = None
        self.acsf = None
        projects = s.query(db.Experiment.project_name).distinct().all()
        project_list = [{'name': str(record[0]), 'type': 'bool'} for record in projects]
        acsf = s.query(db.Experiment.acsf).distinct().all()
        acsf_list = [{'name': str(record[0]), 'type': 'bool'} for record in acsf]
        internal = s.query(db.Experiment.internal).distinct().all()
        internal_list = [{'name': str(record[0]), 'type': 'bool'} for record in internal]
        self.params = Parameter.create(name='Data Filters', type='group', children=[
            {'name': 'Projects', 'type': 'group', 'children':project_list},
            {'name': 'ACSF', 'type': 'group', 'children':acsf_list, 'expanded': False},
            {'name': 'Internal', 'type': 'group', 'children': internal_list, 'expanded': False},
        ])
        self.params.sigTreeStateChanged.connect(self.invalidate_output)

    def get_pair_list(self, session):
        """ Given a set of user selected experiment filters, return a list of pairs.
        Internally uses multipatch_analysis.connectivity.query_pairs.
        """
        if self.pairs is None:
            project_names = [child.name() for child in self.params.child('Projects').children() if child.value() is True]
            project_names = project_names if len(project_names) > 0 else None 
            acsf_recipes = [child.name() for child in self.params.child('ACSF').children() if child.value() is True]
            acsf_recipes = acsf_recipes if len(acsf_recipes) > 0 else None 
            internal_recipes = [child.name() for child in self.params.child('Internal').children() if child.value() is True]
            internal_recipes = internal_recipes if len(internal_recipes) > 0 else None 
            self.pairs = query_pairs(project_name=project_names, acsf=acsf_recipes, session=session, internal=internal_recipes).all()
        return self.pairs

    def invalidate_output(self):
        self.pairs = None
        self.sigOutputChanged.emit(self)

class CellClassFilter(object):
    def __init__(self, cell_class_groups):
        self.cell_groups = None
        self.cell_classes = None
        self._signalHandler = SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged
        self.experiment_filter = ExperimentFilter() 
        self.cell_class_groups = cell_class_groups.keys()
        cell_group_list = [{'name': group, 'type': 'bool'} for group in self.cell_class_groups]
        self.params = Parameter.create(name="Cell Classes", type="group", children=cell_group_list)

        self.params.sigTreeStateChanged.connect(self.invalidate_output)

    def get_cell_groups(self, pairs):
        """Given a list of cell pairs, return a dict indicating which cells
        are members of each user selected cell class.
        This internally calls cell_class.classify_cells
        """
        if self.cell_groups is None:
            self.cell_classes = []
            for group in self.params.children():
                if group.value() is True:
                    self.cell_classes.extend(cell_class_groups[group.name()]) 
            self.cell_classes = [CellClass(**c) for c in self.cell_classes]
            self.cell_groups = classify_cells(self.cell_classes, pairs=pairs)
        return self.cell_groups, self.cell_classes

    def invalidate_output(self):
        self.cell_groups = None
        self.cell_classes = None

class MatrixDisplayFilter(object):
    def __init__(self, view_box, data_fields, confidence_fields):
        self.output = None
        self._signalHandler = SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged
        self.view_box = view_box
        self.legend = None
        self.colorMap = ColorMapParameter()

        self.colorMap.setFields(data_fields)

        self.params = Parameter.create(name='Matrix Display', type='group', children=[
            self.colorMap,
            {'name': 'Text format', 'type': 'str'},
            {'name': 'Show Confidence', 'type': 'list', 'values': confidence_fields},
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

    def element_display_output(self, result):
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
            mappable_result = {k:v for k,v in result.items() if np.isscalar(v)}
            color = colormap.map(mappable_result)[0]
            
            # desaturate low confidence cells
            if result[show_confidence] is not None:
                lower, upper = result[show_confidence]
                confidence = (1.0 - (upper - lower)) ** 2
                color = color * confidence + default_bgcolor * (1.0 - confidence)
        
            # invert text color for dark background
            self.output['fgcolor'] = 'w' if sum(color[:3]) < 384 else 'k'
            self.output['text'] = text_format.format(**result)
            self.output['bgcolor'] = tuple(color)

        return self.output


    def invalidate_output(self):
        self.output = None

# class DataScatterPlot(object):
#     def __init__(self):
#         self.output = None
#         self._signalHandler = SignalHandler()
#         self.sigOutputChanged = self._signalHandler.sigOutputChanged
#         self.spw = pg.ScatterPlotWidget()

#     def set_fields(self, data_fields):
#         self.spw.setFields(data_fields)

#     def set_plot_data(self, data):
#         self.spw.setData(data)
        
#     def invalidate_output(self):
#         self.output = None

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


class ScatterPlot(pg.GraphicsLayoutWidget):
    def __init__(self):
        pg.GraphicsLayoutWidget.__init__(self)
        self.setRenderHints(self.renderHints() | pg.QtGui.QPainter.Antialiasing)


class TracePlot(pg.GraphicsLayoutWidget):
    def __init__(self):
        pg.GraphicsLayoutWidget.__init__(self)
        self.setRenderHints(self.renderHints() | pg.QtGui.QPainter.Antialiasing)


class MatrixAnalyzer(object):
    def __init__(self, session, default_preset=None):
        self.win = MainWindow()
        self.win.show()
        self.win.setGeometry(280, 130,1500, 900)
        self.win.setWindowTitle('Matrix Analyzer')
        self.analyzers = [ConnectivityAnalyzer(), StrengthAnalyzer(), DynamicsAnalyzer()]
        self.active_analyzers = []
        self.scatter_plot = None
        self.line = None
        self.scatter = None
        self.trace_plot = None
        self.trace_plot_list = []
        self.element = None
        self.selected = 0
        self.colors = [(0, 255, 0), (255, 0, 0), (0, 0, 255), (254, 169, 0), (170, 0, 127), (0, 230, 230)]

        self.output_fields, self.confidence_fields = get_all_output_fields(self.analyzers)
        # self.data_scatter_plot = self.win.dsp
        # self.data_scatter_plot.set_fields(self.output_fields)
        # self.spw = self.data_scatter_plot.spw
        self.experiment_filter = ExperimentFilter()
        self.cell_class_filter = CellClassFilter(cell_class_groups)
        self.matrix_display_filter = MatrixDisplayFilter(self.win.matrix_widget.view_box, self.output_fields, self.confidence_fields)
        self.visualizers = [self.matrix_display_filter]
        self.default_preset = default_preset
        self.session = session
        self.results = {}
        self.cell_groups = None
        self.cell_classes = None
        self.params = ptree.Parameter.create(name='params', type='group', children=[
            self.experiment_filter.params, 
            self.cell_class_filter.params,
            self.matrix_display_filter.params,
        ])
        self.win.ptree.setParameters(self.params, showTop=False)

        self.field_map = {}
        for analyzer in self.analyzers:
            for field in analyzer.output_fields()['color_by']:
                self.field_map[field[0]] = analyzer
            for field in analyzer.output_fields()['show_confidence']:
                if field == 'None':
                    continue
                self.field_map[field] = analyzer

        self.presets = self.analyzer_presets()
        preset_list = [p['name'] for p in self.presets]
        preset_params = Parameter.create(name='Analyzer Presets', type='list', values=preset_list)
        self.params.insertChild(0, preset_params)
        
        self.win.update_button.clicked.connect(self.update_clicked)
        self.win.matrix_widget.sigClicked.connect(self.display_matrix_element_data)
        self.experiment_filter.sigOutputChanged.connect(self.cell_class_filter.invalidate_output)
        self.params.child('Analyzer Presets').sigValueChanged.connect(self.presetChanged)
        # self.matrix_display_filter.sigOutputChanged.connect(self.analyzers_needed) 

        # connect up analyzers
        for analyzer in self.analyzers:
            for visualizer in self.visualizers:
                analyzer.sigOutputChanged.connect(visualizer.invalidate_output)
            self.cell_class_filter.sigOutputChanged.connect(analyzer.invalidate_output)

    def analyzer_presets(self):
        self.presets = [
            {'name': 'None'},
            {'name': 'Mouse full connectivity matrix', 'filters': {
            'Data Filters': {
                'Projects': {
                'mouse V1 coarse matrix': True,
            }},
            'Cell Classes': {
                'Mouse All Cre-types by layer': True,
            },
            'Matrix Display': {
                'color_by': 'connection_probability',
                'Text format': '{n_connected}/{n_probed}',
                'Show Confidence': 'cp_confidence_interval',
                'log_scale': True,
            },
            # 'Scatter Plot': [
            #     'connection_probability',
            # ]
            }},
            {'name': 'testing', 'filters': {
                'Data Filters': {
                    'Projects': {
                    'mouse V1 coarse matrix': True,
                }},
                'Cell Classes': {
                    'Mouse Layer 2/3': True,
                },
                'Matrix Display': {
                    'color_by': 'ic_fit_amp_all',
                    'Text format': '{ic_fit_amp_all.mV}',
                    'Show Confidence': 'None',
                    'log_scale': False,
            },
            }},
            ]

        return self.presets

    def presetChanged(self):
        self.clear_preset_selections()
        selected = self.params['Analyzer Presets']
        self.set_preset_selections(selected)
        
    def clear_preset_selections(self):
        for field in self.experiment_filter.params.children():
            for item in field.children():
                item.setValue(False)

        [item.setValue(False) for item in self.cell_class_filter.params.children()]
        # [item.setSelected(False) for item in self.spw.fieldList.selectedItems()]

        self.matrix_display_filter.params.child('Color Map').clearChildren()
        self.matrix_display_filter.params.child('Text format').setValue('')
        self.matrix_display_filter.params.child('Show Confidence').setValue('None')
        self.matrix_display_filter.params.child('log_scale').setValue(False)

    def set_preset_selections(self, selected):
        if selected != 'None':
            selected_preset = [preset for preset in self.presets if preset['name'] == selected][0] 
            for filt, field in selected_preset['filters'].items():
                if filt == 'Scatter Plot':
                    self.spw.setSelectedFields(*field)
                else:
                    for field_name, value in field.items():
                        if isinstance(value, dict):
                            for subfield_name, v2 in value.items():
                                self.params.child(filt).child(field_name).child(subfield_name).setValue(v2)
                        elif field_name == 'color_by':
                            self.matrix_display_filter.colorMap.addNew(value)
                        else:
                            self.params.child(filt).child(field_name).setValue(value)

    def analyzers_needed(self):
        
        data_needed = set()
        for metric in self.matrix_display_filter.colorMap.children():
            data_needed.add(metric.name())
        text_fields = re.findall('\{(.*?)\}', self.matrix_display_filter.params['Text format'])
        for metric in text_fields:
            if ':' in metric:
                data_needed.add(metric.split(':')[0])
            else:
                data_needed.add(metric)
        data_needed.add(self.matrix_display_filter.params['Show Confidence'])
        
        # for metric in self.spw.selectedItems():
        #     data_needed.add(str(metric.txt()))

        
        analyzers = set([self.field_map.get(field, None) for field in data_needed])
        self.active_analyzers = [analyzer for analyzer in analyzers if analyzer is not None]

        self.data_needed = data_needed

        # analyzer_set = set()
        # output_fields, _ = get_all_output_fields(available_analyzers)
        # output_fields = [field[0] for field in output_fields]
        # for data_field in data_needed:
        #     if data_field in output_fields:
        #         analyzer_set.add(analyzer)

        # remove_analyzers = list(set(self.active_analyzers) - analyzer_set)
        # for analyzer in remove_analyzers:
        #     # analyzer.sigOutputChanged.disconnect(self.matrix_display_filter.invalidate_output)
        #     # self.cell_class_filter.sigOutputChanged.disconnect(analyzer.invalidate_output)
        #     self.active_analyzers.remove(analyzer)

        # add_analyzers = list(analyzer_set - set(self.active_analyzers))
        # for analyzer in add_analyzers:
        #     analyzer.sigOutputChanged.connect(self.matrix_display_filter.invalidate_output)
        #     self.cell_class_filter.sigOutputChanged.connect(analyzer.invalidate_output)
        #     self.active_analyzers.append(analyzer)

        print ('Active analyzers:')
        print (self.active_analyzers)

        return self.active_analyzers

    def update_clicked(self):
        with pg.BusyCursor():
            self.analyzers_needed()
            self.update_matrix_results()
            self.update_matrix_display()
            self.display_element_reset()

    def display_matrix_element_data(self, matrix_widget, event, row, col):
        pre_class, post_class = self.matrix_map[row, col]
        analyzer = self.field_map[self.field_name]
        data = analyzer.print_element_info(pre_class, post_class, self.field_name)
        if self.scatter_plot is not None:
            if int(event.modifiers() & pg.QtCore.Qt.ControlModifier)>0:
                self.selected += 1
                if self.selected >= len(self.colors):
                    self.selected = 0
                self.display_element_output(row, col, data, analyzer, trace_plot_list=self.trace_plot_list)
            else:
                self.display_element_reset() 
                self.display_element_output(row, col, data, analyzer)

    def display_element_output(self, row, col, data, analyzer, trace_plot_list=None):
        color = self.colors[self.selected]
        self.element = self.win.matrix_widget.matrix.cells[row][col]
        self.element.setPen(pg.mkPen({'color': color, 'width': 5}))
        pre_class, post_class = self.matrix_map[row, col]
        self.trace_plot = self.win.trace_plot.addPlot()
        self.trace_plot_list.append(self.trace_plot)
        self.line, self.scatter = analyzer.plot_element_data(pre_class, post_class, self.field_name, data=data, color=color, trace_plt=self.trace_plot)
        if len(self.trace_plot_list) > 1:
            first_plot = self.trace_plot_list[0]
            for plot in self.trace_plot_list[1:]:
                plot.setYLink(first_plot)
        self.scatter_plot.addItem(self.line)
        if self.scatter is not None:
            self.scatter_plot.addItem(self.scatter)

    def display_element_reset(self):
        self.selected = 0
        if self.scatter_plot is not None:
            [self.scatter_plot.removeItem(item) for item in self.scatter_plot.items[1:]]
        if self.trace_plot is not None:
           self.win.trace_plot.clear()
           self.trace_plot = None
        self.update_matrix_display()
        self.trace_plot_list = []

    def update_matrix_results(self):
        # Select pairs (todo: age, acsf, internal, temp, etc.)
        self.pairs = self.experiment_filter.get_pair_list(self.session)

        # Group all cells by selected classes
        self.cell_groups, self.cell_classes = self.cell_class_filter.get_cell_groups(self.pairs)
        

        # Group pairs into (pre_class, post_class) groups
        self.pair_groups = classify_pairs(self.pairs, self.cell_groups)

        # analyze matrix elements
        for i, analysis in enumerate(self.active_analyzers):
            results = analysis.measure(self.pair_groups)
            if i == 0:
                self.results = results
            else:
                for pair_group, data in results.items():
                    self.results[pair_group].update(data)
        for pair_group, data in self.results.items():
            no_data = any([data.get('conn_no_data', None), data.get('strength_no_data', None), data.get('dynamics_no_data', None)])
            self.results[pair_group]['no_data'] = no_data

        # set data in scatter plot widget

    def update_matrix_display(self):
        
        shape = (len(self.cell_groups),) * 2
        text = np.empty(shape, dtype=object)
        fgcolor = np.empty(shape, dtype=object)
        bgcolor = np.empty(shape, dtype=object)
        bordercolor = np.empty(shape, dtype=object)
        self.matrix_display_filter.colormap_legend()

        # call display function on every matrix element
        
        self.matrix_map = {}
        for i,row in enumerate(self.cell_groups):
            for j,col in enumerate(self.cell_groups):
                output = self.matrix_display_filter.element_display_output(self.results[(row, col)])
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
        for i,cell_class in enumerate(self.cell_classes):
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

        self.win.matrix_widget.set_matrix_data(text=text, fgcolor=fgcolor, bgcolor=bgcolor, border_color=bordercolor,
                    rows=rows, cols=rows, size=50, header_color='k')

        # plot hist or scatter of data in side panel, if there are multiple colormaps take the first enabled one
        color_map_fields = self.matrix_display_filter.colorMap.children()
        if len(color_map_fields) > 1:
            self.field_name = [field.name() for field in color_map_fields if field['Enabled'] is True][0]
        else:
            self.field_name = color_map_fields[0].name()
        field = self.matrix_display_filter.colorMap.fields[self.field_name]
        if self.scatter_plot is not None:
            self.win.scatter_plot.removeItem(self.scatter_plot)
        self.scatter_plot = self.win.scatter_plot.addPlot()
        results_scatter(self.results, self.field_name, field, self.scatter_plot)

        # self.analysis.summary(self.results, self.field_name)
        


if __name__ == '__main__':

    import sys
    import pyqtgraph as pg
    app = pg.mkQApp()
    pg.dbg()

    session = db.Session()

    
    # Define cell classes
    cell_class_groups = OrderedDict([
        ('Mouse All Cre-types by layer', [
            {'cre_type': 'unknown', 'target_layer': '2/3'},
            #{'pyramidal': True, 'target_layer': '2/3'},
            {'cre_type': 'pvalb', 'target_layer': '2/3'},
            {'cre_type': 'sst', 'target_layer': '2/3'},
            {'cre_type': 'vip', 'target_layer': '2/3'},
           # {'cre_type': 'rorb', 'target_layer': '4'},
            {'cre_type': 'nr5a1', 'target_layer': '4'},
            {'cre_type': 'pvalb', 'target_layer': '4'},
            {'cre_type': 'sst', 'target_layer': '4'},
            {'cre_type': 'vip', 'target_layer': '4'},
            {'cre_type': 'sim1', 'target_layer': '5'},
            {'cre_type': 'tlx3', 'target_layer': '5'},
            {'cre_type': 'pvalb', 'target_layer': '5'},
            {'cre_type': 'sst', 'target_layer': '5'},
            {'cre_type': 'vip', 'target_layer': '5'},
            {'cre_type': 'ntsr1', 'target_layer': '6'},
            {'cre_type': 'pvalb', 'target_layer': '6'},
            {'cre_type': 'sst', 'target_layer': '6'},
            {'cre_type': 'vip', 'target_layer': '6'},
        ]),

        ('Mouse Layer 2/3', [
            {'cre_type': 'unknown', 'target_layer': '2/3'},
            #{'pyramidal': True, 'target_layer': '2/3'},
            {'cre_type': 'pvalb', 'target_layer': '2/3'},
            {'cre_type': 'sst', 'target_layer': '2/3'},
            {'cre_type': 'vip', 'target_layer': '2/3'},
        ]),
        
        ('Mouse Layer 4', [
            {'cre_type': 'nr5a1', 'target_layer': '4'},
            {'cre_type': 'pvalb', 'target_layer': '4'},
            {'cre_type': 'sst', 'target_layer': '4'},
            {'cre_type': 'vip', 'target_layer': '4'},
        ]),

        ('Mouse Layer 5', [
            {'cre_type': ('sim1', 'fam84b'), 'target_layer': '5', 'display_names': ('L5', 'PT\nsim1, fam84b')},
            {'cre_type': 'tlx3', 'target_layer': '5', 'display_names': ('L5', 'IT\ntlx3')},
            {'cre_type': 'pvalb', 'target_layer': '5'},
            {'cre_type': 'sst', 'target_layer': '5'},
            {'cre_type': 'vip', 'target_layer': '5'},
        ]),

        ('Mouse Layer 6', [
            {'cre_type': 'ntsr1', 'target_layer': '6'},
            {'cre_type': 'pvalb', 'target_layer': '6'},
            {'cre_type': 'sst', 'target_layer': '6'},
            {'cre_type': 'vip', 'target_layer': '6'},
        ]),

        ('Mouse Inhibitory Cre-types',[
            {'cre_type': 'pvalb'},
            {'cre_type': 'sst'},
            {'cre_type': 'vip'},
        ]),
 
        ('Mouse Excitatory Cre-types', [
            # {'pyramidal': True, 'target_layer': '2/3'},
            {'cre_type': 'unknown', 'target_layer': '2/3'},
            {'cre_type': 'nr5a1', 'target_layer': '4'},
            {'cre_type': 'sim1', 'target_layer': '5'},
            {'cre_type': 'tlx3', 'target_layer': '5'},
            {'cre_type': 'ntsr1', 'target_layer': '6'},
        ]),

        ('Mouse E-I Cre-types by layer',[
            # {'pyramidal': True, 'target_layer': '2/3'},
            {'cre_type': 'unknown', 'target_layer': '2/3'},
            {'cre_type': ('pvalb', 'sst', 'vip'), 'target_layer': '2/3', 'display_names': ('L2/3', 'Inhibitory\npvalb, sst, vip')},
            {'cre_type': 'nr5a1', 'target_layer': '4'},
            {'cre_type': ('pvalb', 'sst', 'vip'), 'target_layer': '4', 'display_names': ('L4', 'Inhibitory\npvalb, sst, vip')},
            {'cre_type': 'sim1', 'target_layer': '5'},
            {'cre_type': 'tlx3', 'target_layer': '5'},
            {'cre_type': ('pvalb', 'sst', 'vip'), 'target_layer': '5', 'display_names': ('L5', 'Inhibitory\npvalb, sst, vip')},
            {'cre_type': 'ntsr1', 'target_layer': '6'},
            {'cre_type': ('pvalb', 'sst', 'vip'), 'target_layer': '6', 'display_names': ('L6', 'Inhibitory\npvalb, sst, vip')},     
        ]),

        ('Pyramidal / Nonpyramidal by layer', [
            {'pyramidal': True, 'target_layer': '2'},
            {'pyramidal': False, 'target_layer': '2'},
            {'pyramidal': True, 'target_layer': '3'},
            {'pyramidal': False, 'target_layer': '3'},
            {'pyramidal': True, 'target_layer': '4'},
            {'pyramidal': False, 'target_layer': '4'},
            {'pyramidal': True, 'target_layer': '5'},
            {'pyramidal': False, 'target_layer': '5'},
            {'pyramidal': True, 'target_layer': '6'},
            {'pyramidal': False, 'target_layer': '6'},
        ]),

        ('Pyramidal by layer', [
            {'pyramidal': True, 'target_layer': '2'}, 
            {'pyramidal': True, 'target_layer': '3'},
            {'pyramidal': True, 'target_layer': '4'},
            {'pyramidal': True, 'target_layer': '5'},
            {'pyramidal': True, 'target_layer': '6'},
        ]),

        ('All cells by layer', [
            {'target_layer': '2'},
            {'target_layer': '3'},
            {'target_layer': '4'},
            {'target_layer': '5'},
            {'target_layer': '6'},
        ]),
    ])


    maz = MatrixAnalyzer(session=session, default_preset='None')

    if sys.flags.interactive == 0:
        app.exec_()