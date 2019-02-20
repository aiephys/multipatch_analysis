"""
Prototype code for analyzing connectivity and synaptic properties between cell classes.


"""

from __future__ import print_function, division

from collections import OrderedDict
import numpy as np
import pyqtgraph as pg
from multipatch_analysis.database import database as db
from multipatch_analysis.connectivity import query_pairs, ConnectivityAnalyzer, StrengthAnalyzer
from multipatch_analysis.connection_strength import ConnectionStrength, get_amps, get_baseline_amps
from multipatch_analysis.morphology import Morphology
from multipatch_analysis import constants
from multipatch_analysis.cell_class import CellClass, classify_cells, classify_pairs
from multipatch_analysis.ui.graphics import MatrixItem
from pyqtgraph import parametertree as ptree
from pyqtgraph.parametertree import Parameter


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
        self.matrix_widget = MatrixWidget()
        self.h_splitter.addWidget(self.matrix_widget)
        self.plot_splitter = pg.QtGui.QSplitter()
        self.plot_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.h_splitter.addWidget(self.plot_splitter)
        self.scatter_plot = ScatterPlot()
        self.trace_plot = TracePlot()
        self.plot_splitter.addWidget(self.scatter_plot)
        self.plot_splitter.addWidget(self.trace_plot)
        self.h_splitter.setSizes([150, 300, 150])


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
        self.params = Parameter.create(name='Data Filters', type='group', children=[
            {'name': 'Projects', 'type': 'group', 'children':project_list},
            {'name': 'ACSF', 'type': 'group', 'children':acsf_list},
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
            self.pairs = query_pairs(project_name=project_names, acsf=acsf_recipes, session=session).all()
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


class Analysis(object):
    def __init__(self):
        self.analyzer = None
        self._signalHandler = SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged
        self.analyzers = {'Connectivity': ConnectivityAnalyzer, 'Strength': StrengthAnalyzer}
        analyzer_list = [{'name': 'Analysis', 'type': 'list', 'values': self.analyzers.keys()}]
        self.params = Parameter.create(name='Analyzers', type='group', children=analyzer_list)

        self.params.sigTreeStateChanged.connect(self.invalidate_output)

    def get_analyzer(self):
        if self.analyzer is None:
            analyzer_key = self.params.children()[0].value()
            self.analyzer = self.analyzers[analyzer_key]
        return self.analyzer

    def invalidate_output(self):
        self.analyzer = None


class DisplayFilter(object):
    def __init__(self):
        self.displays = None
        self._signalHandler = SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged
        self.params = Parameter.create(name='Display Options', type='group', children=[
            {'name': 'Matrix Color Gradient', 'type': 'colormap', 'value': pg.ColorMap(
            [0, 0.01, 0.03, 0.1, 0.3, 1.0],
            [(0,0,100), (80,0,80), (140,0,0), (255,100,0), (255,255,100), (255,255,255)])},
            {'name': 'Shade by CI', 'type': 'bool', 'value': True},
            {'name': 'Periodic Table Style', 'type': 'bool', 'value': False},
        ])

        self.params.sigTreeStateChanged.connect(self.invalidate_output)
        
        # def get_displays(self):
        #     if self.displays is None:
        #         self.disp

    def invalidate_output(self):
        self.displays = None


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
    def __init__(self, session):
        self.session = session

        self.experiment_filter = ExperimentFilter()
        self.cell_class_filter = CellClassFilter(cell_class_groups)
        self.analysis = Analysis()
        self.display_filter = DisplayFilter()
        self.params = ptree.Parameter.create(name='params', type='group', children=[
            self.experiment_filter.params, 
            self.cell_class_filter.params, 
            self.analysis.params, 
            self.display_filter.params
        ])

        self.win = MainWindow()
        self.win.show()
        self.win.resize(1500, 800)
        self.win.setWindowTitle('Matrix Analyzer')
        self.win.ptree.setParameters(self.params, showTop=False)

        self.win.update_button.clicked.connect(self.update_clicked)
        self.win.matrix_widget.sigClicked.connect(self.display_matrix_element_data)
        self.experiment_filter.sigOutputChanged.connect(self.cell_class_filter.invalidate_output)
        self.cell_class_filter.sigOutputChanged.connect(self.analysis.invalidate_output)
        self.analysis.sigOutputChanged.connect(self.display_filter.invalidate_output)

    def update_clicked(self):
        with pg.BusyCursor():
            self.update_matrix()

    def element_connection_list(self, pre_class, post_class):
        results = ConnectivityAnalyzer(self.pair_groups).measure()
        connections = results[(pre_class, post_class)]['connected_pairs']
        print ("Connection type: %s -> %s" % (pre_class, post_class))
        print ("Connected Pairs:")
        for connection in connections:
            print ("\t %s" % (connection))
        probed_pairs = results[(pre_class, post_class)]['probed_pairs']
        print ("Probed Pairs:")
        for probed in probed_pairs:
            print ("\t %s" % (probed))

    def display_matrix_element_data(self, matrix_widget, event, row, col):
        pre_class, post_class = self.matrix_map[row, col]
        self.element_connection_list(pre_class, post_class)

    def update_matrix(self):
        # Select pairs (todo: age, acsf, internal, temp, etc.)
        self.pairs = self.experiment_filter.get_pair_list(self.session)

        # Group all cells by selected classes
        cell_groups, cell_classes = self.cell_class_filter.get_cell_groups(self.pairs)
        

        # Group pairs into (pre_class, post_class) groups
        self.pair_groups = classify_pairs(self.pairs, cell_groups)

        # analyze matrix elements
        analyzer = self.analysis.get_analyzer()
        matrix_analysis = analyzer(self.pair_groups)
        results = matrix_analysis.measure()

        shape = (len(cell_groups),) * 2
        text = np.empty(shape, dtype=object)
        fgcolor = np.empty(shape, dtype=object)
        bgcolor = np.empty(shape, dtype=object)
        bordercolor = np.empty(shape, dtype=object)

        # call display function on every matrix element
        
        self.matrix_map = {}
        for i,row in enumerate(cell_groups):
            for j,col in enumerate(cell_groups):
                output = matrix_analysis.display(row, col, results[(row, col)])
                self.matrix_map[i, j] = (row, col)
                text[i, j] = output['text']
                fgcolor[i, j] = output['fgcolor']
                bgcolor[i, j] = output['bgcolor']
                bordercolor[i, j] = output['bordercolor']
                
        # Force cell class descriptions down to tuples of 2 items
        # Kludgy, but works for now.
        rows = []
        for cell_class in cell_classes:
            tup = cell_class.as_tuple
            row = tup[:1]
            if len(tup) > 1:
                row = row + (' '.join(tup[1:]),)
            rows.append(row)

        self.win.matrix_widget.set_matrix_data(text=text, fgcolor=fgcolor, bgcolor=bgcolor, border_color=bordercolor,
                    rows=rows, cols=rows, size=50, header_color='k')
        
        matrix_analysis.summary(results)


if __name__ == '__main__':

    import pyqtgraph as pg
    pg.dbg()

    session = db.Session()
    
    # Define cell classes
    cell_class_groups = {
        'Mouse All Cre-types by layer': [
            # {'cre_type': 'unknown', 'pyramidal': True, 'target_layer': '2/3'},
            # {'cre_type': 'unknown', 'target_layer': '2/3'},
            # {'pyramidal': True, 'target_layer': '2/3'},
            {'pyramidal': True, 'target_layer': '2/3'},
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
        ],

        'Pyramidal / Nonpyramidal by layer': [
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
        ],
    }


    maz = MatrixAnalyzer(session=session)
