import multipatch_analysis.database as db
import pyqtgraph as pg
import sys
from pyqtgraph import parametertree as ptree
from pyqtgraph.parametertree import Parameter
from pyqtgraph.widgets.DataFilterWidget import DataFilterParameter
from neuroanalysis.ui.plot_grid import PlotGrid
from multipatch_analysis.ui.experiment_browser import ExperimentBrowser
from collections import OrderedDict
from neuroanalysis.data import Trace, TraceList

default_latency = 11e-3

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
        # self.next_pair_button = pg.QtGui.QPushButton("Load Next Pair")
        # self.v_splitter.addWidget(self.next_pair_button)
        self.h_splitter.addWidget(self.vc_plot.grid)
        self.h_splitter.addWidget(self.ic_plot.grid)
        self.h_splitter.setSizes([200, 175, 175])
        self.layout.addWidget(self.h_splitter)
        
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
        self.do_fit = Parameter.create(name='Fit PSP', type='action')
        self.fit_params = Parameter.create(name='Fit parameters', type='group', children=[
            {'name': 'Parameter', 'type': 'str', 'value': '-55 VC, -70 VC, -55 IC, -70 IC', 'readonly': True},
            {'name': 'Amplitude', 'type': 'float', 'readonly': True},
            {'name': 'Latency', 'type': 'float', 'readonly': True},
            {'name': 'Rise time', 'type': 'float', 'readonly': True},
            {'name': 'Decay tau', 'type': 'float', 'readonly': True},
            ])
        self.fit_pass = Parameter.create(name="Fit Pass", type='bool')
        self.save_params = Parameter.create(name='Save Analysis', type='action')
        self.comments = Parameter.create(name='Comments', type='text')
        self.params = Parameter.create(name='params', type='group', children=[
            self.latency,
            self.synapse,
            self.gap,
            self.do_fit,
            self.fit_params,
            self.fit_pass,
            self.save_params,
            self.comments
            ])
        
    def update_params(self, **kargs):
        for k, v in kargs.items():
            self.params.child(k).setValue(v)

    def set_ic_latency(self, ic_superline):
        value = ic_superline.pos()
        self.params.child('Latency', 'IC').setValue(value)

    def set_vc_latency(self, vc_superline):
        value = vc_superline.pos()
        self.params.child('Latency', 'VC').setValue(value)

class VCPlot(pg.GraphicsLayoutWidget):
    def __init__(self, superline):
        pg.GraphicsLayoutWidget.__init__(self)
        self.grid = PlotGrid()
        self.grid.set_shape(2, 1)
        self.grid.show()
        self.plots = (self.grid[0, 0], self.grid[1, 0])
        self.plots[0].setTitle('Voltage Clamp')
        self.plots[0].hideAxis('bottom')
        self.plots[0].addItem(superline.new_line(default_latency))
        self.plots[1].addItem(superline.new_line(default_latency))
        self.plots[1].setXRange(5e-3, 20e-3)
        self.plots[0].setXLink(self.plots[1])

    def plot_traces(self, traces_dict):
        for i, traces in enumerate(traces_dict.values()):
            if len(traces) == 0:
                continue
            for trace in traces:
                self.plots[i].plot(trace.time_values, trace.data, pen=(255,255,255,10))
            grand_trace = TraceList(traces).mean()
            self.plots[i].plot(grand_trace.time_values, grand_trace.data, pen={'color': 'b', 'width': 2})

    def clear_plots(self):
        self.plots[0].clear()
        self.plots[1].clear()

class ICPlot(pg.GraphicsLayoutWidget):
    def __init__(self, superline):
        pg.GraphicsLayoutWidget.__init__(self)
        self.grid = PlotGrid()
        self.grid.set_shape(2, 1)
        self.grid.show()
        self.plots = (self.grid[0, 0], self.grid[1, 0])
        self.plots[0].setTitle('Current Clamp')
        self.plots[0].hideAxis('bottom')
        self.plots[0].addItem(superline.new_line(default_latency))
        self.plots[1].addItem(superline.new_line(default_latency))
        self.plots[1].setXRange(5e-3, 30e-3)
        self.plots[0].setXLink(self.plots[1])

    def plot_traces(self, traces_dict):  
        for i, traces in enumerate(traces_dict.values()):
            if len(traces) == 0:
                continue
            for trace in traces:
                self.plots[i].plot(trace.time_values, trace.data, pen=(255,255,255,100))
            grand_trace = TraceList(traces).mean()
            self.plots[i].plot(grand_trace.time_values, grand_trace.data, pen={'color': 'b', 'width': 2})

    def clear_plots(self):
        self.plots[0].clear()
        self.plots[1].clear()

class SuperLine(pg.QtCore.QObject):
    sigPositionChanged = pg.QtCore.Signal(object)

    def __init__(self):
        pg.QtCore.QObject.__init__(self)
        self.lines = []

    def new_line(self, x_pos):
        self.line = pg.InfiniteLine(x_pos, pen={'color': 'y', 'width': 3}, movable=True)
        self.line.setZValue(100)
        self.lines.append(self.line)
        self.line.sigPositionChanged.connect(self.line_sync)
        return self.line

    def line_sync(self, moved_line):
        for line in self.lines:
            with pg.SignalBlock(line.sigPositionChanged, self.line_sync):
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
        self.ic_plot = ICPlot(self.ic_superline)
        self.vc_superline = SuperLine()
        self.vc_superline.sigPositionChanged.connect(self.ctrl_panel.set_vc_latency)
        self.vc_plot = VCPlot(self.vc_superline)
        self.experiment_browser = ExperimentBrowser()
            # start with latency of 11 ms, fit, then user can adjust if needed
            # connect fit button to fitting algorithm (which one? from where?)
                # need to pull set latency value as starting point
            # plot fit on top of traces
            # connect save button to dump result into new DB
            # connect load next pair button to clear results and move to the next pair

    def load_pair(self, pair):
        print ('loading responses...')
        s = db.Session()
        q = self.response_query(s, pair)
        self.pulse_responses = q.all()
        print('got this many responses: %d' % len(self.pulse_responses))
        s.close()

        pair_params = {'Synapse call': pair.connection_strength.synapse_type, 'Gap junction call': pair.electrical}
        self.ctrl_panel.update_params(**pair_params)

        self.analyze_responses()


    def analyze_responses(self):
        ex_limits = [-80e-3, -63e-3]
        in_limits = [-62e-3, -45e-3]
        self.grand_trace = OrderedDict([('vc', {'-55': [], '-70': []}), ('ic', {'-55': [], '-70': []})])
        for rec in self.pulse_responses:
            if rec.ind_freq not in [10, 20, 50]:
                continue
            data = rec.data
            start_time = rec.rec_start
            spike_time = rec.spike_time
            clamp = rec.clamp_mode
            holding = rec.baseline_potential
            data_trace = Trace(data=data, t0=start_time-spike_time+10e-3, sample_rate=db.default_sample_rate)
            if clamp == 'vc':
                if in_limits[0] < holding < in_limits[1]:
                    self.grand_trace['vc']['-55'].append(data_trace)
                elif ex_limits[0] < holding < ex_limits[1]:
                    self.grand_trace['vc']['-70'].append(data_trace)
            if clamp == 'ic':
                if in_limits[0] < holding < in_limits[1]:
                    self.grand_trace['ic']['-55'].append(data_trace)
                elif ex_limits[0] < holding < ex_limits[1]:
                    self.grand_trace['ic']['-70'].append(data_trace)        
            
        self.vc_plot.plot_traces(self.grand_trace['vc'])
        self.ic_plot.plot_traces(self.grand_trace['ic'])


    def response_query(self, session, pair):
        q = session.query(
        db.PulseResponse.id.label('response_id'),
        db.PulseResponse.data,
        db.PulseResponse.start_time.label('rec_start'),
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


if __name__ == '__main__':

    app = pg.mkQApp()
    pg.dbg()

    s = db.Session()
    q = s.query(db.Experiment).filter(db.Experiment.acq_timestamp==1511913615.871)
    expts = q.all()

    mw = MainWindow()
    mw.set_expts(expts)


    if sys.flags.interactive == 0:
        app.exec_()