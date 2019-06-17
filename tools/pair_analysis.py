import multipatch_analysis.database as db
import pyqtgraph as pg
import sys
import numpy as np
from pyqtgraph import parametertree as ptree
from pyqtgraph.parametertree import Parameter
from pyqtgraph.widgets.DataFilterWidget import DataFilterParameter
from neuroanalysis.ui.plot_grid import PlotGrid
from multipatch_analysis.ui.experiment_browser import ExperimentBrowser
from collections import OrderedDict
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.baseline import float_mode
from multipatch_analysis.connection_detection import fit_psp
from random import shuffle

default_latency = 11e-3
comment_hashtag = [
    '',
    '#doublespike',
    '#badspikes',
    '#fixable',
    '#secondopinion',
    '#lostcause']

class SignalHandler(pg.QtCore.QObject):
        """Because we can't subclass from both QObject and QGraphicsRectItem at the same time
        """
        sigOutputChanged = pg.QtCore.Signal(object) #self

class MainWindow(pg.QtGui.QWidget):
    def __init__(self):
        pg.QtGui.QWidget.__init__(self)
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
        self.ic_plot = self.pair_analyzer.ic_plot
        self.vc_plot = self.pair_analyzer.vc_plot
        self.experiment_browser = self.pair_analyzer.experiment_browser
        self.v_splitter = pg.QtGui.QSplitter()
        self.v_splitter.setOrientation(pg.QtCore.Qt.Vertical)
        self.h_splitter.addWidget(self.v_splitter)
        self.v_splitter.addWidget(self.experiment_browser)
        self.v_splitter.addWidget(self.ptree)
        self.v_splitter.setSizes([100, 300])
        # self.next_pair_button = pg.QtGui.QPushButton("Load Next Pair")
        # self.v_splitter.addWidget(self.next_pair_button)
        self.h_splitter.addWidget(self.vc_plot.grid)
        self.h_splitter.addWidget(self.ic_plot.grid)
        self.h_splitter.setSizes([200, 200, 200])
        self.layout.addWidget(self.h_splitter)
        self.setGeometry(280, 130, 1500, 900)
        self.show()

        # self.next_pair_button.clicked.connect(self.load_next_pair)
        self.experiment_browser.itemSelectionChanged.connect(self.selected_pair)
        

    def set_expts(self, expts):
        self.experiment_browser.populate(experiments=expts)

    def selected_pair(self):
        selected = self.experiment_browser.selectedItems()
        if len(selected) != 1:
            return
        item = selected[0]
        if hasattr(item, 'pair') is False:
            return
        pair = item.pair
        self.pair_param.setValue(pair)
        self.pair_analyzer.load_pair(pair)


class ControlPanel(object):
    def __init__(self):
        self.latency = Parameter.create(name='Latency', type='group', children=[
            {'name': 'VC', 'type': 'float', 'suffix': 's', 'siPrefix': True, 'readonly': True, 'value': default_latency},
            {'name': 'IC', 'type': 'float', 'suffix': 's', 'siPrefix': True, 'readonly': True, 'value': default_latency},
            ])
        self.synapse = Parameter.create(name='Synapse call', type='list', values={'Excitatory': 'ex', 'Inhibitory': 'in', 'None': None})
        self.gap = Parameter.create(name='Gap junction call', type='bool')
        fit_cat = ['-55 VC', '-70 VC', '-55 IC', '-70 IC']
        fit_param = [{'name': c, 'type': 'group', 'children': [
            {'name': 'Parameter:', 'type': 'str', 'readonly': True, 'value': 'Amplitude, Latency, Rise time, Decay tau, NRMSE'},
            {'name': 'Value:', 'type': 'str', 'readonly': True},
            {'name': 'Fit Pass', 'type': 'bool'},
            ]} for c in fit_cat]
        self.fit_params = Parameter.create(name='Fit parameters', type='group', children=fit_param)
        self.save_params = Parameter.create(name='Save Analysis', type='action')
        self.comments = Parameter.create(name='Comments', type='group', children=[
            {'name': 'Hashtag', 'type': 'list', 'values': comment_hashtag},
            {'name': '', 'type': 'text'}])
        self.params = Parameter.create(name='params', type='group', children=[
            self.latency,
            self.synapse,
            self.gap,
            self.fit_params,
            self.save_params,
            self.comments
            ])

        self.comments.child('Hashtag').sigValueChanged.connect(self.add_text_to_comments)

    def add_text_to_comments(self):
        text = self.comments['Hashtag']
        comments = self.comments[''] + '/n'
        update_comments = comments + text
        self.comments.child('').setValue(update_comments)
        
    def update_params(self, **kargs):
        for k, v in kargs.items():
            self.params.child(k).setValue(v)

    def update_fit_params(self, fit_params):
        param_names = ['amp', 'xoffset', 'rise_time', 'decay_tau', 'nrmse']
        for mode in ['vc', 'ic']:
            if mode == 'vc':
                suffix = ['A', 's', 's', 's', '']
            elif mode == 'ic':
                suffix = ['V', 's', 's', 's', '']
            for holding in ['-55', '-70']:
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

class TracePlot(pg.GraphicsLayoutWidget):
    def __init__(self):
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
        self.qc_color = {'qc_pass': (255,255,255,100), 'qc_fail': (150, 0, 0)}

    def plot_traces(self, traces_dict):  
        for i, holding in enumerate(traces_dict.keys()):
            for qc, traces in traces_dict[holding].items():
                if len(traces) == 0:
                    continue
                for trace in traces:
                    self.trace_plots[i].plot(trace.time_values, trace.data, pen=self.qc_color[qc])
                if qc == 'qc_pass':
                    grand_trace = TraceList(traces).mean()
                    self.trace_plots[i].plot(grand_trace.time_values, grand_trace.data, pen={'color': 'b', 'width': 2})
            self.trace_plots[i].autoRange()
            self.trace_plots[i].setXRange(5e-3, 20e-3)
            # y_range = [grand_trace.data.min(), grand_trace.data.max()]
            # self.plots[i].setYRange(y_range[0], y_range[1], padding=1)

    def plot_spikes(self, spikes_dict):
        for i, holding in enumerate(spikes_dict.keys()):
            for qc, spikes in spikes_dict[holding].items():
                if len(spikes) == 0:
                    continue
                for spike in spikes:
                    self.spike_plots[i].plot(spike.time_values, spike.data, pen=self.qc_color[qc])

    def plot_fit(self, trace, fit, holding, fit_pass):
        if holding == '-55':
            self.trace_plots[0].addLegend()
            if self.fit_item_55 is not None:
                self.trace_plots[0].removeItem(self.fit_item_55)
            self.fit_item_55 = pg.PlotDataItem(trace.time_values, fit.best_fit, name='-55 holding', pen={'color': self.fit_color[False], 'width': 3})
            self.trace_plots[0].addItem(self.fit_item_55)
        
        elif holding == '-70':
            self.trace_plots[1].addLegend()  
            if self.fit_item_70 is not None:
                self.trace_plots[1].removeItem(self.fit_item_70)
            self.fit_item_70 = pg.PlotDataItem(trace.time_values, fit.best_fit, name='-70 holding', pen={'color': self.fit_color[False], 'width': 3})
            self.trace_plots[1].addItem(self.fit_item_70)

    def color_fit(self, name, value):
        if '-55' in name:
            if self.fit_item_55 is not None:
                self.fit_item_55.setPen({'color': self.fit_color[value], 'width': 3})
        if '-70' in name:
            if self.fit_item_70 is not None:
                self.fit_item_70.setPen({'color': self.fit_color[value], 'width': 3})

    def clear_plots(self):
        for plot in self.plots:
            plot.clear()
        
        self.plots[-1].autoRange()
        self.plots[-1].setXRange(5e-3, 20e-3)
        self.fit_item_70 = None
        self.fit_item_55 = None


class VCPlot(TracePlot):
    def __init__(self, superline):
        TracePlot.__init__(self)
        self.plots[0].setTitle('Voltage Clamp')
        for plot in self.plots[:-1]:
            plot.addItem(superline.new_line(default_latency))
            plot.setXLink(self.plots[-1])
        self.plots[-1].setXRange(5e-3, 20e-3)
        self.plots[-1].addItem(superline.new_line(default_latency))

        for plot in self.plots:
            plot.setLabel('left', units='A')

class ICPlot(TracePlot):
    def __init__(self, superline):
        TracePlot.__init__(self)
        self.plots[0].setTitle('Current Clamp')
        for plot in self.plots[:-1]:
            plot.addItem(superline.new_line(default_latency))
            plot.setXLink(self.plots[-1])
        self.plots[-1].setXRange(5e-3, 20e-3)
        self.plots[-1].addItem(superline.new_line(default_latency))

        for plot in self.plots:
            plot.setLabel('left', units='V')
        

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

class PairAnalysis(object):
    def __init__(self):
        self.ctrl_panel = ControlPanel()
        self.ic_superline = SuperLine()
        self.ic_superline.sigPositionChanged.connect(self.ctrl_panel.set_ic_latency)
        self.ic_superline.sigPositionChangeFinished.connect(self.fit_response_update)
        self.ic_plot = ICPlot(self.ic_superline)
        self.vc_superline = SuperLine()
        self.vc_superline.sigPositionChanged.connect(self.ctrl_panel.set_vc_latency)
        self.vc_superline.sigPositionChangeFinished.connect(self.fit_response_update)
        self.vc_plot = VCPlot(self.vc_superline)
        self.ctrl_panel.params.child('Fit parameters').sigTreeStateChanged.connect(self.colorize_fit)
        self.experiment_browser = ExperimentBrowser()
        self.fit_params = {}
        self.fit_pass = False
        self.nrmse_thresh = 4
        self.use_x_range = True
        self.signs = {'vc': {
        '-55': {'ex': '-', 'in': '+'}, 
        '-70': {'ex': '-', 'in': 'any'},
        },
        'ic':{
        '-55': {'ex': '+', 'in': '-'},
        '-70': {'ex': '+', 'in': 'any'},
        },
        }
            # connect save button to dump result into new DB
            # connect load next pair button to clear results and move to the next pair

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
        for p, plot in enumerate(self.vc_plot.plots):
            line = self.vc_superline.lines[p]
            line.setValue(default_latency)
            plot.addItem(line)
        for p, plot in enumerate(self.ic_plot.plots):
            line = self.ic_superline.lines[p]
            line.setValue(default_latency)
            plot.addItem(line)
        self.ctrl_panel.params.child('Comments', 'Hashtag').setValue('')
        self.ctrl_panel.params.child('Comments', '').setValue('')
        
    def load_pair(self, pair):
        with pg.BusyCursor():
            self.pair = pair
            self.reset_display()
            print ('loading responses...')
            s = db.Session()
            q = self.response_query(s, pair)
            self.pulse_responses = q.all()
            print('got this many responses: %d' % len(self.pulse_responses))
            s.close()

            if pair.synapse is True:
                synapse_type = pair.connection_strength.synapse_type
            else:
                synapse_type = None
            pair_params = {'Synapse call': synapse_type, 'Gap junction call': pair.electrical}
            self.ctrl_panel.update_params(**pair_params)

            self.analyze_responses()

    def analyze_responses(self):
        ex_limits = [-80e-3, -63e-3]
        in_limits = [-62e-3, -45e-3]
        qc = {False: 'qc_fail', True: 'qc_pass'}
        self.traces = OrderedDict([('vc', {'-55': {'qc_pass': [], 'qc_fail': []}, '-70': {'qc_pass': [], 'qc_fail': []}}), 
                                ('ic', {'-55': {'qc_pass': [], 'qc_fail': []}, '-70': {'qc_pass': [], 'qc_fail': []}})])
        self.spikes = OrderedDict([('vc', {'-55': {'qc_pass': [], 'qc_fail': []}, '-70': {'qc_pass': [], 'qc_fail': []}}), 
                                ('ic', {'-55': {'qc_pass': [], 'qc_fail': []}, '-70': {'qc_pass': [], 'qc_fail': []}})])
        for rec in self.pulse_responses:
            if rec.ind_freq not in [10, 20, 50]:
                continue
            data = rec.data
            spike = rec.spike
            n_spikes = rec.n_spikes
            start_time = rec.rec_start
            spike_time = rec.spike_time if rec.spike_time is not None else 0. 
            clamp = rec.clamp_mode
            holding = rec.baseline_potential
            t0 = start_time-spike_time+10e-3
            baseline = float_mode(data[0:int(db.default_sample_rate*6e-3)])
            data_trace = Trace(data=data-baseline, t0=t0, sample_rate=db.default_sample_rate)
            spike_trace = Trace(data=spike, t0=t0, sample_rate=db.default_sample_rate)
            trace_qc_pass = rec.ex_qc_pass
            spike_qc_pass = n_spikes == 1
            
            if in_limits[0] < holding < in_limits[1]:
                self.traces[clamp]['-55'][qc[trace_qc_pass]].append(data_trace)
                self.spikes[clamp]['-55'][qc[spike_qc_pass]].append(spike_trace)
            elif ex_limits[0] < holding < ex_limits[1]:
                self.traces[clamp]['-70'][qc[trace_qc_pass]].append(data_trace)
                self.spikes[clamp]['-70'][qc[spike_qc_pass]].append(spike_trace)

        self.vc_plot.plot_traces(self.traces['vc'])
        self.vc_plot.plot_spikes(self.spikes['vc'])
        self.ic_plot.plot_traces(self.traces['ic'])
        self.ic_plot.plot_spikes(self.spikes['ic'])

        self.fit_responses()

    def response_query(self, session, pair):
        q = session.query(
        db.PulseResponse.id.label('response_id'),
        db.PulseResponse.data,
        db.PulseResponse.ex_qc_pass,
        db.PulseResponse.start_time.label('rec_start'),
        db.StimPulse.data.label('spike'),
        db.StimPulse.n_spikes,
        db.StimSpike.max_dvdt_time.label('spike_time'),
        db.PatchClampRecording.clamp_mode,
        db.PatchClampRecording.baseline_potential,
        db.MultiPatchProbe.induction_frequency.label('ind_freq'),
        )
        q = q.join(db.StimPulse, db.PulseResponse.stim_pulse)
        q = q.join(db.StimSpike, db.StimSpike.stim_pulse_id==db.StimPulse.id)
        q = q.join(db.Recording, db.PulseResponse.recording)
        q = q.join(db.PatchClampRecording)
        q = q.join(db.MultiPatchProbe)
        q = q.filter(db.PulseResponse.pair_id==pair.id)

        return q

    def fit_response_update(self):
        self.use_x_range = False
        self.fit_responses(x_offset_win=[-0.1e-3, 0.1e-3])

    def fit_responses(self, x_offset_win=[-1e-3, 6e-3]):
        initial_fit_parameters = OrderedDict([('vc', {'-55': {}, '-70': {}}), ('ic', {'-55': {}, '-70': {}})])
        output_fit_parameters = OrderedDict([('vc', {'-55': {}, '-70': {}}), ('ic', {'-55': {}, '-70': {}})])
        for mode in ['vc', 'ic']:
            for holding in ['-55', '-70']:
                self.fit_pass = False
                self.ctrl_panel.params.child('Fit parameters', holding + ' ' + mode.upper(), 'Fit Pass').setValue(self.fit_pass)
                sign = self.signs[mode][holding].get(self.ctrl_panel.params['Synapse call'], 'any')
                x_offset = self.ctrl_panel.params['Latency', mode.upper()]
                initial_fit_parameters[mode][holding]['xoffset'] = x_offset
                if len(self.traces[mode][holding]['qc_pass']) == 0:
                    continue
                grand_trace = TraceList(self.traces[mode][holding]['qc_pass']).mean()
                base_rgn = grand_trace.time_slice(-6e-3, 0)
                weight = np.ones(len(grand_trace.data))*10.  #set everything to ten initially
                weight[int(12e-3/db.default_sample_rate):int(19e-3/db.default_sample_rate)] = 30.  #area around steep PSP rise 
                if mode == 'vc':
                    stacked = False
                    initial_rise = 1e-3
                    rise_bounds = [0.1e-3, 5e-3]
                elif mode == 'ic':
                    stacked  = True
                    initial_rise = 5e-3
                    rise_bounds = [1e-3, 25e-3]
                    weight[int(10e-3/db.default_sample_rate):int(12e-3/db.default_sample_rate)] = 0.   #area around stim artifact
                x_win = [x_offset + x_offset_win[0], x_offset + x_offset_win[1]]
                rise_times = list(initial_rise*2.**np.arange(-2, 3, 0.5))
                if self.use_x_range is True:
                    x_range = (list(np.linspace(x_win[0], x_win[1], 7)), x_win[0], x_win[1])
                else:
                    x_range = (x_offset, 'fixed')
                try:
                    fit = fit_psp(grand_trace, 
                        mode=mode, 
                        sign=sign, 
                        xoffset=x_range, 
                        rise_time=(rise_times, rise_bounds[0], rise_bounds[1]),
                        stacked=stacked,
                        )
                    for param, val in fit.best_values.items():
                        if param == 'xoffset':
                            val  = val - 10e-3
                        output_fit_parameters[mode][holding][param] = val
                    output_fit_parameters[mode][holding]['yoffset'] = fit.best_values['yoffset']
                    output_fit_parameters[mode][holding]['nrmse'] = fit.nrmse()
                    self.fit_pass = fit.nrmse() < self.nrmse_thresh
                    if mode == 'vc':
                        self.vc_plot.plot_fit(grand_trace, fit, holding, self.fit_pass)
                    elif mode == 'ic':
                        self.ic_plot.plot_fit(grand_trace, fit, holding, self.fit_pass)
                    self.ctrl_panel.params.child('Fit parameters', holding + ' ' + mode.upper(), 'Fit Pass').setValue(self.fit_pass)
                except:
                    print("Error in PSP fit:")
                    sys.excepthook(*sys.exc_info())
                    continue
                
        self.fit_params = {'initial': initial_fit_parameters, 'fit': output_fit_parameters}
        self.ctrl_panel.update_fit_params(output_fit_parameters)


if __name__ == '__main__':

    app = pg.mkQApp()
    pg.dbg()

    # expt_list = [1532552839.296, 1536184462.404, 1534802671.146, 1532552954.427, 1557178633.726, 1517266234.424, 1554243759.900]
    s = db.Session()
    e = s.query(db.Experiment.acq_timestamp)
    expt_list = e.all()
    expt_list = [ee[0] for ee in expt_list]
    shuffle(expt_list)
    expt_list = expt_list[:5]
    q = s.query(db.Experiment).filter(db.Experiment.acq_timestamp.in_(expt_list))
    expts = q.all()

    mw = MainWindow()
    mw.set_expts(expts)


    if sys.flags.interactive == 0:
        app.exec_()