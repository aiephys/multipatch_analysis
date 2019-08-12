from multipatch_analysis.database import default_db as db
import multipatch_analysis.data_notes_db as notes_db
import pyqtgraph as pg
import sys, copy, argparse
import numpy as np
from pyqtgraph import parametertree as ptree
from pyqtgraph.parametertree import Parameter
from pyqtgraph.widgets.DataFilterWidget import DataFilterParameter
from neuroanalysis.ui.plot_grid import PlotGrid
from multipatch_analysis.ui.experiment_browser import ExperimentBrowser
from collections import OrderedDict
from neuroanalysis.data import TSeriesList
from neuroanalysis.fitting import Psp, StackedPsp
from multipatch_analysis.avg_response_fit import response_query, sort_responses, fit_avg_response, pair_notes_query
from random import shuffle, seed

default_latency = 1e-3
comment_hashtag = [
    '#doublespike',
    '#doublepsp',
    '#badspikes',
    '#fixable',
    '#needsecondopinion',
    '#lostcause',
    '#MVP',
    '#crosstalk',
    '#badqc',
    '#risetime',
    '#baseline',
    '#nodatatofit',
    '#polysynaptic',
    '#doublechecked',
    '#syncretypmismatch',
]

comment_hashtag.sort(key=lambda x:x[1])
comment_hashtag = [''] + comment_hashtag

modes = ['vc', 'ic']
holdings = ['-55', '-70']


class MainWindow(pg.QtGui.QWidget):
    def __init__(self, default_session, notes_session):
        pg.QtGui.QWidget.__init__(self)
        self.default_session = default_session
        self.notes_session = notes_session
        self.layout = pg.QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.h_splitter = pg.QtGui.QSplitter()
        self.h_splitter.setOrientation(pg.QtCore.Qt.Horizontal)
        self.pair_analyzer = PairAnalysis()
        self.ctrl_panel = self.pair_analyzer.ctrl_panel
        self.params = self.ctrl_panel.params
        self.ptree = ptree.ParameterTree(showHeader=False)
        self.pair_param = Parameter.create(name='Current Pair', type='str', readonly=True)
        self.ptree.addParameters(self.pair_param)
        self.ptree.addParameters(self.params, showTop=False)
        self.save_btn = pg.FeedbackButton('Save Analysis')
        self.expt_btn = pg.QtGui.QPushButton('Set Experiments with Hashtags')
        self.ic_plot = self.pair_analyzer.ic_plot
        self.vc_plot = self.pair_analyzer.vc_plot
        self.hash_ptree = ptree.ParameterTree(showHeader=False)
        self.hash_select = Parameter.create(name='Select Hashtags', type='group', children=
            [{'name': 'With multiple selected:', 'type': 'list', 'values': ['Include if any appear', 'Include if all appear'], 'value': 'Include if any appear'}]+
            [{'name': '#', 'type': 'bool'}] +
            [{'name': ht, 'type': 'bool'} for ht in comment_hashtag[1:]])
        self.hash_ptree.addParameters(self.hash_select)
        self.experiment_browser = self.pair_analyzer.experiment_browser
        self.v_splitter = pg.QtGui.QSplitter()
        self.v_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.h_splitter.addWidget(self.v_splitter)
        self.v_splitter.addWidget(self.hash_ptree)
        self.v_splitter.addWidget(self.expt_btn)
        self.v_splitter.addWidget(self.experiment_browser)
        self.v_splitter.addWidget(self.ptree)
        self.v_splitter.addWidget(self.save_btn)
        self.v_splitter.setSizes([50, 20, 100, 300, 20])
        # self.next_pair_button = pg.QtGui.QPushButton("Load Next Pair")
        # self.v_splitter.addWidget(self.next_pair_button)
        self.h_splitter.addWidget(self.vc_plot.grid)
        self.h_splitter.addWidget(self.ic_plot.grid)
        self.fit_compare = self.pair_analyzer.fit_compare
        self.meta_compare = self.pair_analyzer.meta_compare
        self.v2_splitter = pg.QtGui.QSplitter()
        self.v2_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.v2_splitter.addWidget(self.fit_compare)
        self.v2_splitter.addWidget(self.meta_compare)
        self.h_splitter.addWidget(self.v2_splitter)
        self.h_splitter.setSizes([200, 200, 200, 500])
        self.layout.addWidget(self.h_splitter)
        self.fit_compare.hide()
        self.meta_compare.hide()
        self.setGeometry(280, 130, 1500, 900)
        self.show()

        # self.next_pair_button.clicked.connect(self.load_next_pair)
        self.experiment_browser.itemSelectionChanged.connect(self.selected_pair)
        self.save_btn.clicked.connect(self.save_to_db)
        self.expt_btn.clicked.connect(self.get_expts_hashtag)

    def save_to_db(self):
        try:
            self.pair_analyzer.save_to_db()
            self.save_btn.success()
        except:
            self.save_btn.failure()
            raise

    def get_expts_hashtag(self):
        selected_hashtags = [ht.name() for ht in self.hash_select.children()[1:] if ht.value() is True]
        q = self.notes_session.query(notes_db.PairNotes)
        pairs_to_include = []
        note_pairs = q.all()
        note_pairs.sort(key=lambda p: p.expt_id)
        for p in note_pairs:
            comments = p.notes['comments']    
            if len(selected_hashtags) == 1:
                hashtag = selected_hashtags[0]
                if hashtag == '#':
                    if hashtag in comments and all([ht not in comments for ht in comment_hashtag[1:]]):
                        print(p.expt_id, p.pre_cell_id, p.post_cell_id, comments)
                        pairs_to_include.append(p)
                else:
                    if hashtag in comments:
                        print(p.expt_id, p.pre_cell_id, p.post_cell_id, comments)
                        pairs_to_include.append(p)
                
            if len(selected_hashtags) > 1:
                hashtag_present = [ht in comments for ht in selected_hashtags]
                or_expts = self.hash_select['With multiple selected:'] == 'Include if any appear'
                and_expts = self.hash_select['With multiple selected:'] == 'Include if all appear'
                if or_expts and any(hashtag_present):
                    print(p.expt_id, p.pre_cell_id, p.post_cell_id, comments)
                    pairs_to_include.append(p)
                if and_expts and all(hashtag_present):
                    print(p.expt_id, p.pre_cell_id, p.post_cell_id, comments)
                    pairs_to_include.append(p)

        timestamps = set([pair.expt_id for pair in pairs_to_include])
        q2 = self.default_session.query(db.Experiment).filter(db.Experiment.acq_timestamp.in_(timestamps))
        expts = q2.all()
        self.set_expts(expts)

    def set_expts(self, expts):
        self.experiment_browser.clear()
        self.experiment_browser.populate(experiments=expts)

    def selected_pair(self):
        self.fit_compare.hide()
        self.meta_compare.hide()
        selected = self.experiment_browser.selectedItems()
        if len(selected) != 1:
            return
        item = selected[0]
        if hasattr(item, 'pair') is False:
            return
        pair = item.pair
        ## check to see if the pair has already been analyzed
        expt_id = '%0.3f' % pair.experiment.acq_timestamp
        pre_cell_id = str(pair.pre_cell.ext_id)
        post_cell_id = str(pair.post_cell.ext_id)
        q = pair_notes_query(self.notes_session, pair)
        record = q.all()
        self.pair_param.setValue(pair)
        if len(record) == 0:
            self.pair_analyzer.load_pair(pair, self.default_session)
            self.pair_analyzer.analyze_responses()
            self.pair_analyzer.fit_responses()
        elif len(record) == 1:
            msg = pg.QtGui.QMessageBox.question(self, "Pair Analysis", 
                "Pair %s %s->%s has already been analyzed. \n Would you like to load the results?" % (expt_id, pre_cell_id, post_cell_id),
                pg.QtGui.QMessageBox.Yes | pg.QtGui.QMessageBox.No)
            if msg == pg.QtGui.QMessageBox.Yes:
                self.pair_analyzer.load_pair(pair, self.default_session, record=record[0])
                self.pair_analyzer.analyze_responses()
                self.pair_analyzer.load_saved_fit(record[0])
        else:
            raise Exception('More than one record for this pair %s %s->%s was found in the Pair Notes database' % (expt_id, pre_cell_id, post_cell_id))


class ControlPanel(object):
    def __init__(self):
        self.latency = Parameter.create(name='Latency', type='group', children=[
            {'name': 'VC', 'type': 'float', 'suffix': 's', 'siPrefix': True, 'readonly': True, 'value': default_latency},
            {'name': 'IC', 'type': 'float', 'suffix': 's', 'siPrefix': True, 'readonly': True, 'value': default_latency},
        ])
        self.synapse = Parameter.create(name='Synapse call', type='list', values={'Excitatory': 'ex', 'Inhibitory': 'in', 'None': None})
        self.gap = Parameter.create(name='Gap junction call', type='bool')
        fit_cat = ['-55 VC', '-70 VC', '-55 IC', '-70 IC']
        fit_param = [
            {'name': c, 'type': 'group', 'children': [
                {'name': 'Parameter:', 'type': 'str', 'readonly': True, 'value': 'Amplitude, Latency, Rise time, Decay tau, NRMSE'},
                {'name': 'Value:', 'type': 'str', 'readonly': True},
                {'name': 'Fit Pass', 'type': 'bool'},
            ]} for c in fit_cat
        ]
        self.fit_params = Parameter.create(name='Fit parameters', type='group', children=fit_param)
        self.warn_param = Parameter.create(name='Warnings', type='text', readonly=True)
        self.comments = Parameter.create(name='Comments', type='group', children=[
            {'name': 'Hashtag', 'type': 'list', 'values': comment_hashtag, 'value': ''},
            {'name': '', 'type': 'text'}
        ])
        self.params = Parameter.create(name='params', type='group', children=[
            self.latency,
            self.synapse,
            self.gap,
            self.fit_params,
            self.warn_param,
            self.comments
        ])

        self.comments.child('Hashtag').sigValueChanged.connect(self.add_text_to_comments)

    def add_text_to_comments(self):
        text = self.comments['Hashtag']
        comments = self.comments['']
        update_comments = comments + text + '\n'
        self.comments.child('').setValue(update_comments)
        
    def update_params(self, **kargs):
        for k, v in kargs.items():
            self.params.child(k).setValue(v)

    def update_fit_params(self, fit_params):
        param_names = ['amp', 'xoffset', 'rise_time', 'decay_tau', 'nrmse']
        for mode in modes:
            if mode == 'vc':
                suffix = ['A', 's', 's', 's', '']
            elif mode == 'ic':
                suffix = ['V', 's', 's', 's', '']
            for holding in holdings:
                group = self.params.child('Fit parameters', holding + ' ' + mode.upper())
                values = fit_params[mode][holding]
                format_values = self.format_fit_output(values, param_names, suffix)
                group.child('Value:').setValue(format_values)

    def format_fit_output(self, values, name, suffix):
        format_list = []
        for p in zip(name, suffix):
            if values.get(p[0]) is None:
                format_list.append('nan')
            else:
                value = values[p[0]]
                if p[0] == 'nrmse':
                    p_format = ('%0.2f' % value)
                else:
                    p_format = pg.siFormat(value, suffix=p[1])
                format_list.append(p_format)
        output = ", ".join(format_list)
        return output

    def set_ic_latency(self, ic_superline):
        value = ic_superline.pos()
        self.params.child('Latency', 'IC').setValue(value)

    def set_vc_latency(self, vc_superline):
        value = vc_superline.pos()
        self.params.child('Latency', 'VC').setValue(value)


class TSeriesPlot(pg.GraphicsLayoutWidget):
    def __init__(self, title, units):
        pg.GraphicsLayoutWidget.__init__(self)
        self.grid = PlotGrid()
        self.grid.set_shape(4, 1)
        self.grid.grid.ci.layout.setRowStretchFactor(0, 3)
        self.grid.grid.ci.layout.setRowStretchFactor(1, 8)
        self.grid.grid.ci.layout.setRowStretchFactor(2, 5)
        self.grid.grid.ci.layout.setRowStretchFactor(3, 10)
        self.grid.show()
        self.trace_plots = (self.grid[1, 0], self.grid[3, 0])
        self.spike_plots = (self.grid[0, 0], self.grid[2, 0])
        self.plots = self.spike_plots + self.trace_plots
        for plot in self.plots[:-1]:
            plot.hideAxis('bottom')
        self.plots[-1].setLabel('bottom', text='Time from spike', units='s')
        self.fit_item_55 = None
        self.fit_item_70 = None
        self.fit_color = {True: 'g', False: 'r'}
        self.qc_color = {'qc_pass': (255, 255, 255, 100), 'qc_fail': (255, 0, 0, 100)}

        self.plots[0].setTitle(title)
        for (plot, holding) in zip(self.trace_plots, holdings):
            plot.setXLink(self.plots[-1])
            plot.setLabel('left', text="%d holding" % int(holding), units=units)
        for plot in self.spike_plots:
            plot.setXLink(self.plots[-1])
            plot.setLabel('left', text="presynaptic spike")
            plot.addLine(x=0)
        self.plots[-1].setXRange(-5e-3, 10e-3)
        
        self.items = []

    def plot_traces(self, traces_dict):  
        for i, holding in enumerate(traces_dict.keys()):
            for qc, traces in traces_dict[holding].items():
                if len(traces) == 0:
                    continue
                for trace in traces:
                    item = self.trace_plots[i].plot(trace.time_values, trace.data, pen=self.qc_color[qc])
                    if qc == 'qc_fail':
                        item.setZValue(-10)
                    self.items.append(item)
                if qc == 'qc_pass':
                    grand_trace = TSeriesList(traces).mean()
                    item = self.trace_plots[i].plot(grand_trace.time_values, grand_trace.data, pen={'color': 'b', 'width': 2})
                    self.items.append(item)
            self.trace_plots[i].autoRange()
            self.trace_plots[i].setXRange(-5e-3, 10e-3)
            # y_range = [grand_trace.data.min(), grand_trace.data.max()]
            # self.plots[i].setYRange(y_range[0], y_range[1], padding=1)

    def plot_spikes(self, spikes_dict):
        for i, holding in enumerate(spikes_dict.keys()):
            for qc, spikes in spikes_dict[holding].items():
                if len(spikes) == 0:
                    continue
                for spike in spikes:
                    item = self.spike_plots[i].plot(spike.time_values, spike.data, pen=self.qc_color[qc])
                    if qc == 'qc_fail':
                        item.setZValue(-10)
                    self.items.append(item)

    def plot_fit(self, trace, fit, holding, fit_pass=False):
        if holding == '-55':
            if self.fit_item_55 is not None:
                self.trace_plots[0].removeItem(self.fit_item_55)
            self.fit_item_55 = pg.PlotDataItem(trace.time_values, fit, name='-55 holding', pen={'color': self.fit_color[fit_pass], 'width': 3})
            self.trace_plots[0].addItem(self.fit_item_55)
        
        elif holding == '-70':
            if self.fit_item_70 is not None:
                self.trace_plots[1].removeItem(self.fit_item_70)
            self.fit_item_70 = pg.PlotDataItem(trace.time_values, fit, name='-70 holding', pen={'color': self.fit_color[fit_pass], 'width': 3})
            self.trace_plots[1].addItem(self.fit_item_70)

    def color_fit(self, name, value):
        if '-55' in name:
            if self.fit_item_55 is not None:
                self.fit_item_55.setPen({'color': self.fit_color[value], 'width': 3})
        if '-70' in name:
            if self.fit_item_70 is not None:
                self.fit_item_70.setPen({'color': self.fit_color[value], 'width': 3})

    def clear_plots(self):
        for item in self.items + [self.fit_item_55, self.fit_item_70]:
            if item is None:
                continue
            item.scene().removeItem(item)
        self.items = []
        
        self.plots[-1].autoRange()
        self.plots[-1].setXRange(-5e-3, 10e-3)
        self.fit_item_70 = None
        self.fit_item_55 = None
        

class SuperLine(pg.QtCore.QObject):
    sigPositionChanged = pg.QtCore.Signal(object)
    sigPositionChangeFinished = pg.QtCore.Signal(object)

    def __init__(self):
        pg.QtCore.QObject.__init__(self)
        self.lines = []

    def new_line(self, x_pos):
        self.line = pg.InfiniteLine(x_pos, pen={'color': 'y', 'width': 3}, movable=True)
        self.line.setZValue(100)
        self.lines.append(self.line)
        self.line.sigPositionChanged.connect(self.line_sync)
        self.line.sigPositionChangeFinished.connect(self.move_finished)
        return self.line

    def move_finished(self):
        self.sigPositionChangeFinished.emit(self)

    def line_sync(self, moved_line):
        for line in self.lines:
            with pg.SignalBlock(line.sigPositionChanged, self.line_sync):
                with pg.SignalBlock(line.sigPositionChangeFinished, self.move_finished):
                    line.setValue(moved_line.value())
        self.sigPositionChanged.emit(self)

    def pos(self):
        line_positions = []
        for line in self.lines:
            line_positions.append(line.value())
        n_positions = set(line_positions)
        if len(n_positions) > 1:
            raise Exception("Lines are out of sync and reporting different positions")
        position = list(n_positions)[0]

        return position
        
    def set_value(self, value, block_fit=False):
        for line in self.lines:
            with pg.SignalBlock(line.sigPositionChanged, self.line_sync):
                with pg.SignalBlock(line.sigPositionChangeFinished, self.move_finished):
                    line.setValue(value)
        self.sigPositionChanged.emit(self)
        if block_fit is False:
            self.sigPositionChangeFinished.emit(self)


class PairAnalysis(object):
    def __init__(self):
        self.ctrl_panel = ControlPanel()
        
        self.ic_superline = SuperLine()
        self.ic_superline.sigPositionChanged.connect(self.ctrl_panel.set_ic_latency)
        self.ic_superline.sigPositionChangeFinished.connect(self.ic_fit_response_update)
        self.ic_plot = TSeriesPlot('Current Clamp', 'V')
        for plot in self.ic_plot.trace_plots:
            plot.addItem(self.ic_superline.new_line(default_latency))

        self.vc_superline = SuperLine()
        self.vc_superline.sigPositionChanged.connect(self.ctrl_panel.set_vc_latency)
        self.vc_superline.sigPositionChangeFinished.connect(self.vc_fit_response_update)
        self.vc_plot = TSeriesPlot('Voltage Clamp', 'A')
        for plot in self.vc_plot.trace_plots:
            plot.addItem(self.vc_superline.new_line(default_latency))
            
        self.ctrl_panel.params.child('Fit parameters').sigTreeStateChanged.connect(self.colorize_fit)
        self.experiment_browser = ExperimentBrowser()
        self.fit_compare = pg.DiffTreeWidget()
        self.meta_compare = pg.DiffTreeWidget()
        self.nrmse_thresh = 4
        self.traces = OrderedDict()
        self.spikes = OrderedDict()
        self.signs = {
            'vc': {
                '-55': {'ex': '-', 'in': '+'}, 
                '-70': {'ex': '-', 'in': 'any'},
            },
            'ic':{
                '-55': {'ex': '+', 'in': '-'},
                '-70': {'ex': '+', 'in': 'any'},
            },
        }

        self.fit_precision = {
            'amp': {'vc': 14, 'ic': 8},
            'exp_amp': {'vc': 14, 'ic': 8},
            'decay_tau': {'vc': 8, 'ic': 8},
            'nrmse': {'vc': 2, 'ic': 2},
            'rise_time': {'vc': 7, 'ic': 6},
            'rise_power': {'vc': 0, 'ic': 0},
            'xoffset': {'vc': 7, 'ic': 7},
            'yoffset': {'vc': 5, 'ic': 14},
        }

    def colorize_fit(self, param, changes):
        for c in changes:
            p, change, info = c
            if p.name() != 'Fit Pass':
                continue
            if 'VC' in p.parent().name():
                self.vc_plot.color_fit(p.parent().name(), info)
            if 'IC' in p.parent().name():
                self.ic_plot.color_fit(p.parent().name(), info)

    def reset_display(self):
        self.vc_plot.clear_plots()
        self.ic_plot.clear_plots()
        self.vc_superline.set_value(default_latency, block_fit=True)
        self.ic_superline.set_value(default_latency, block_fit=True)
        self.ctrl_panel.params.child('Comments', 'Hashtag').setValue('')
        self.ctrl_panel.params.child('Comments', '').setValue('')
        self.ctrl_panel.params.child('Warnings').setValue('')
        
    def load_pair(self, pair, default_session, record=None):
        self.record = hash(record)
        self.initial_fit_parameters = OrderedDict([
            ('vc', {'-55': {}, '-70': {}}), 
            ('ic', {'-55': {}, '-70': {}}),
        ])
        self.output_fit_parameters = OrderedDict([
            ('vc', {'-55': {}, '-70': {}}), 
            ('ic', {'-55': {}, '-70': {}}),
        ])
        self.fit_params = {'initial': self.initial_fit_parameters, 'fit': self.output_fit_parameters}

        with pg.BusyCursor():
            self.reset_display()
            self.pair = pair
            print ('loading responses...')
            q = response_query(default_session, pair)
            self.pulse_responses = q.all()
            print('got this many responses: %d' % len(self.pulse_responses))
                
            if pair.synapse is True:
                synapse_type = pair.connection_strength.synapse_type if pair.connection_strength is not None else None
            else:
                synapse_type = None
            pair_params = {'Synapse call': synapse_type, 'Gap junction call': pair.electrical}
            self.ctrl_panel.update_params(**pair_params)

    def analyze_responses(self):
        self.traces, self.spikes = sort_responses(self.pulse_responses)
        fitable_responses = []
        for mode in modes:
            for holding in holdings:
                fitable_responses.append(bool(self.traces[mode][holding]['qc_pass']))
        if not any(fitable_responses):
            print('No fitable responses, bailing out')
        self.vc_plot.plot_traces(self.traces['vc'])
        self.vc_plot.plot_spikes(self.spikes['vc'])
        self.ic_plot.plot_traces(self.traces['ic'])
        self.ic_plot.plot_spikes(self.spikes['ic'])

    def ic_fit_response_update(self):
        latency = self.ctrl_panel.params['Latency', 'IC']
        self.fit_responses(mode='ic', latency=latency)

    def vc_fit_response_update(self):
        latency = self.ctrl_panel.params['Latency', 'VC']
        self.fit_responses(mode='vc', latency=latency)

    def fit_responses(self, mode=None, latency=None):
        if mode is not None:
            c_modes = [mode]
        else:
            c_modes = modes
        with pg.ProgressDialog("curve fitting..", maximum=len(modes)*len(holdings)) as dlg:
            for mode in c_modes:
                for holding in holdings:
                    self.fit_pass = False
                    sign = self.signs[mode][holding].get(self.ctrl_panel.params['Synapse call'], 'any')
                    ofp, x_offset, best_fit = fit_avg_response(self.traces, mode, holding, latency, sign)
                    self.initial_fit_parameters[mode][holding]['xoffset'] = x_offset
                    self.output_fit_parameters[mode][holding].update(ofp)
                    self.fit_pass = ofp.get('nrmse', self.nrmse_thresh) < self.nrmse_thresh
                    self.ctrl_panel.params.child('Fit parameters', holding + ' ' + mode.upper(), 'Fit Pass').setValue(self.fit_pass)
                    if mode == 'vc'and bool(ofp):
                        avg_trace = TSeriesList(self.traces['vc'][holding]['qc_pass']).mean()
                        self.vc_plot.plot_fit(avg_trace, best_fit, holding, self.fit_pass)
                    elif mode == 'ic' and bool(ofp):
                        avg_trace = TSeriesList(self.traces['ic'][holding]['qc_pass']).mean()
                        self.ic_plot.plot_fit(avg_trace, best_fit, holding, self.fit_pass)
                    dlg += 1
                    if dlg.wasCanceled():
                        raise Exception("User canceled fit")
        self.fit_params['initial'].update(self.initial_fit_parameters)
        self.fit_params['fit'].update(self.output_fit_parameters)
        self.ctrl_panel.update_fit_params(self.fit_params['fit'])
        self.generate_warnings() 

    def generate_warnings(self):
        self.warnings = []
        latency_mode = []
        for mode in modes:
            latency_holding = []
            for holding in holdings:
                initial_latency = self.fit_params['initial'].get(mode).get(holding).get('xoffset')
                fit_latency = self.fit_params['fit'].get(mode).get(holding).get('xoffset')
                if fit_latency is None or initial_latency is None:
                    continue
                if abs(np.diff([fit_latency, initial_latency])) > 0.1e-3:
                    warning = 'Initial latency and fit latency differ by %s' % pg.siFormat(abs(np.diff([fit_latency, initial_latency])), suffix='s')
                latency_holding.append(fit_latency)
            if len(latency_holding) == 2 and abs(np.diff(latency_holding)) > 0.01e-3:
                warning = 'Fit latencies for %s mode do not match' % mode
                self.warnings.append(warning)
            latency_mode.append(np.mean(latency_holding))
        latency_diff = np.diff(latency_mode)[0]
        if  abs(latency_diff) > 0.2e-3:
            warning = 'Latency across modes differs by %s' % pg.siFormat(latency_diff, suffix='s')
            self.warnings.append(warning)

        if np.min(latency_mode) < 0.4e-3 and self.ctrl_panel.params['Gap junction call'] is False:
            self.warnings.append("Short latency; is this a gap junction?")

        guess_sign = []
        guess = None
        for mode, fits1 in self.output_fit_parameters.items():
            for holding, fit in fits1.items():
                if 'amp' not in fit:
                    continue
                if mode == 'ic':
                    guess = 1 if fit['amp'] > 0 else -1
                elif mode == 'vc':
                    guess = -1 if fit['amp'] > 0 else 1
                guess_sign.append(guess)
        if np.all(np.array(guess) == 1):
            guess = "ex"
        elif np.all(np.array(guess) == -1):
            guess = "in"
        if guess is None:
            self.warnings.append("Mixed amplitude signs; pick ex/in carefully.")
        elif guess != self.ctrl_panel.params['Synapse call']:
            self.warnings.append("Looks like an %s synapse??" % guess)

        print_warning = '\n'.join(self.warnings)
        self.ctrl_panel.params.child('Warnings').setValue(print_warning)

    def save_to_db(self):
        fit_pass = {}
        for mode in modes:
            fit_pass[mode] = {}
            for holding in holdings:
                fit_pass[mode][holding] = self.ctrl_panel.params['Fit parameters', holding + ' ' + mode.upper(), 'Fit Pass']

        expt_id = '%0.3f' % self.pair.experiment.acq_timestamp
        pre_cell_id = str(self.pair.pre_cell.ext_id)   
        post_cell_id = str(self.pair.post_cell.ext_id)
        meta = {
            'expt_id': expt_id,
            'pre_cell_id': pre_cell_id,
            'post_cell_id': post_cell_id,
            'synapse_type': self.ctrl_panel.params['Synapse call'],
            'gap_junction': self.ctrl_panel.params['Gap junction call'],
            'fit_parameters': self.fit_params,
            'fit_pass': fit_pass,
            'fit_warnings': self.warnings,
            'comments': self.ctrl_panel.params['Comments', ''],
        }
        
        fields = {
            'expt_id': expt_id,
            'pre_cell_id': pre_cell_id,
            'post_cell_id': post_cell_id, 
            'notes': meta,
        }
        s = notes_db.db.session(readonly=False)
        q = pair_notes_query(s, self.pair)
        rec = q.all()
        if len(rec) == 0:
            record_check = hash(None)
        elif len(rec) == 1:
            saved_rec = rec[0]
            record_check = hash(saved_rec)
        else:
            raise Exception('More than one record was found for pair %s %s->%s in the Pair Notes database' % (expt_id, pre_cell_id, post_cell_id))

        if hash(self.record) == record_check:
            entry = notes_db.PairNotes(**fields)
            s.add(entry)
            s.commit()
        else:
            self.print_pair_notes(meta, saved_rec)
            msg = pg.QtGui.QMessageBox.question(None, "Pair Analysis", 
                "The record you are about to save conflicts with what is in the Pair Notes database.\nYou can see the differences highlighted in red.\nWould you like to overwrite?",
                pg.QtGui.QMessageBox.Yes | pg.QtGui.QMessageBox.No)
            if msg == pg.QtGui.QMessageBox.Yes:
                saved_rec.expt_id = expt_id
                saved_rec.pre_cell_id = pre_cell_id
                saved_rec.post_cell_id = post_cell_id
                saved_rec.notes = meta
                s.commit() 
            else:
                raise Exception('Save Cancelled')
        s.close()

    def print_pair_notes(self, meta, saved_rec):
        meta_copy = copy.deepcopy(meta)
        current_fit = {k:v for k, v in meta_copy['fit_parameters']['fit'].items()}
        saved_fit = {k:v for k, v in saved_rec.notes['fit_parameters']['fit'].items()}
        
        for mode in modes:
            for holding in holdings:
                current_fit[mode][holding] = {k:round(v, self.fit_precision[k][mode]) for k, v in current_fit[mode][holding].items()}
                saved_fit[mode][holding] = {k:round(v, self.fit_precision[k][mode]) for k, v in saved_fit[mode][holding].items()}

        self.fit_compare.setData(current_fit, saved_fit)
        self.fit_compare.trees[0].setHeaderLabels(['Current Fit Parameters', 'type', 'value'])
        self.fit_compare.trees[1].setHeaderLabels(['Saved Fit Parameters', 'type', 'value'])

        current_meta = {k:v for k, v in meta_copy.items() if k != 'fit_parameters'} 
        saved_meta = {k:v for k, v in saved_rec.notes.items() if k != 'fit_parameters'} 
        saved_meta.update({k:str(saved_meta[k]) for k in ['comments', 'expt_id', 'pre_cell_id', 'post_cell_id']})

        self.meta_compare.setData(current_meta, saved_meta)
        self.meta_compare.trees[0].setHeaderLabels(['Current Metadata', 'type', 'value'])
        self.meta_compare.trees[1].setHeaderLabels(['Saved Metadata', 'type', 'value'])

        self.fit_compare.show()
        self.meta_compare.show()

    def load_saved_fit(self, record):
        data = record.notes
        pair_params = {'Synapse call': data['synapse_type'], 'Gap junction call': data['gap_junction']}
        self.ctrl_panel.update_params(**pair_params)
        self.ctrl_panel.update_fit_params(data['fit_parameters']['fit'])
        self.warnings = data['fit_warnings']
        self.ctrl_panel.params.child('Warnings').setValue('\n'.join(self.warnings))
        self.ctrl_panel.params.child('Comments', '').setValue(data['comments'])
        self.fit_params = data['fit_parameters']
        self.output_fit_parameters = data['fit_parameters']['fit']
        self.initial_fit_parameters = data['fit_parameters']['initial']

        initial_vc_latency = data['fit_parameters']['initial']['vc']['-55']['xoffset']
        self.vc_superline.set_value(initial_vc_latency, block_fit=True)
        initial_ic_latency = data['fit_parameters']['initial']['ic']['-55']['xoffset']
        self.ic_superline.set_value(initial_ic_latency, block_fit=True)
        for mode in modes:
            for holding in holdings:
                fit_pass = data['fit_pass'][mode][holding]
                self.ctrl_panel.params.child('Fit parameters', holding + ' ' + mode.upper(), 'Fit Pass').setValue(fit_pass)
                fit_params = copy.deepcopy(data['fit_parameters']['fit'][mode][holding])
                
                if fit_params:
                    fit_params.pop('nrmse', None)
                    if mode == 'ic':
                        p = StackedPsp() 
                    if mode == 'vc':
                        p = Psp()
                        
                    avg = TSeriesList(self.traces[mode][holding]['qc_pass']).mean()
                    
                    fit_psp = p.eval(x=avg.time_values, **fit_params)
                    
                    if mode == 'vc':
                        self.vc_plot.plot_fit(avg, fit_psp, holding, fit_pass=fit_pass)
                    if mode == 'ic':
                        self.ic_plot.plot_fit(avg, fit_psp, holding, fit_pass=fit_pass)            



if __name__ == '__main__':
    from sqlalchemy import or_

    app = pg.mkQApp()
    pg.dbg()
    parser = argparse.ArgumentParser()
    parser.add_argument('--user', type=int)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--hashtag', action='store_true')
    parser.add_argument('--timestamps', type=float, nargs='*')

    args = parser.parse_args(sys.argv[1:])
    user = args.user
    n_users = 10

    default_session = db.session()
    notes_session = notes_db.db.session()
    synapses = default_session.query(db.Pair).filter(db.Pair.synapse==True).all()
    timestamps = set([pair.experiment.acq_timestamp for pair in synapses])
    
    mw = MainWindow(default_session, notes_session)
    if user is not None:
        user_nums = [(ts, int(ts*1000) % n_users) for ts in timestamps]
        timestamps = [un[0] for un in user_nums if un[1] == args.user]
    elif args.check is True:
        pair_in_notes = []
        pair_not_in_notes = []
        ghost_pair = []
        print('checking %d pairs...' % len(synapses))
        for pair in synapses:
            pair_notes = pair_notes_query(notes_session, pair).all()
            if len(pair_notes) == 0:
                pair_not_in_notes.append(pair)
            elif len(pair_notes) == 1:
                pair_in_notes.append(pair)
            else:
                ghost_pair.append(pair)
        timestamps = set([pair.experiment.acq_timestamp for pair in pair_not_in_notes])
        print('%d pairs in notes db' % len(pair_in_notes))
        print('%d pairs not in notes db' % len(pair_not_in_notes))
        print('%d pairs mysteriously missing' % (len(ghost_pair)))
        print('%d/%d pairs accounted for' % (sum([len(pair_in_notes), len(pair_not_in_notes), len(ghost_pair)]), len(synapses)))
    elif args.timestamps is not None:
        timestamps = args.timestamps   
    else:
        seed(0)
        timestamps = list(timestamps)
        shuffle(timestamps)
        timestamps = timestamps[:10]
    
    q = default_session.query(db.Experiment).filter(db.Experiment.acq_timestamp.in_(timestamps))
    expts = q.all()
    mw.set_expts(expts)


    if sys.flags.interactive == 0:
        app.exec_()
