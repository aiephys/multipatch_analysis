import multipatch_analysis.database as db
import pyqtgraph as pg
import sys
from pyqtgraph import parametertree as ptree
from pyqtgraph.parametertree import Parameter
from pyqtgraph.widgets.DataFilterWidget import DataFilterParameter
from neuroanalysis.ui.plot_grid import PlotGrid
from multipatch_analysis.ui.experiment_browser import ExperimentBrowser
from multipatch_analysis.pulse_response_strength import response_query

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
        self.latency = Parameter.create(name='Latency', type='float', readonly=True, suffix='m', siPrefix=True)
        self.synapse = Parameter.create(name='Synapse call', type='list', values={'Excitatory': 'ex', 'Inhibitory': 'in', 'None': None})
        self.gap = Parameter.create(name='Gap junction call', type='bool')
        self.do_fit = Parameter.create(name='Fit PSP', type='action')
        self.fit_params = Parameter.create(name='Fit parameters', type='group', readonly=True, children=[
            {'name': 'Parameter', 'type': 'str', 'value': '-55 VC, -70 VC, -55 IC, -70 IC', 'readonly': True},
            {'name': 'Amplitude', 'type': 'float', 'readonly': True},
            {'name': 'Latency', 'type': 'float', 'readonly': True},
            {'name': 'Rise time', 'type': 'float', 'readonly': True},
            {'name': 'Decay tau', 'type': 'float', 'readonly': True},
            ])
        self.fit_pass = Parameter.create(name="Fit Pass", type='bool')
        self.save_params = Parameter.create(name='Save Analysis', type='action')
        self.params = Parameter.create(name='params', type='group', children=[
            self.latency,
            self.synapse,
            self.gap,
            self.do_fit,
            self.fit_params,
            self.fit_pass,
            self.save_params,
            ])
        

class VCPlot(pg.GraphicsLayoutWidget):
    def __init__(self):
        pg.GraphicsLayoutWidget.__init__(self)
        self.grid = PlotGrid()
        self.grid.set_shape(2, 1)
        self.grid.show()
        self.plots = (self.grid[0, 0], self.grid[1, 0])
        self.plots[0].setTitle('Voltage Clamp')

    def clear_plots(self):
        self.plots[0].clear()
        self.plots[1].clear()

class ICPlot(pg.GraphicsLayoutWidget):
    def __init__(self):
        pg.GraphicsLayoutWidget.__init__(self)
        self.grid = PlotGrid()
        self.grid.set_shape(2, 1)
        self.grid.show()
        self.plots = (self.grid[0, 0], self.grid[1, 0])
        self.plots[0].setTitle('Current Clamp')

    def clear_plots(self):
        self.plots[0].clear()
        self.plots[1].clear()

class PairAnalysis(object):
    def __init__(self):
        self.ctrl_panel = ControlPanel()
        self.ic_plot = ICPlot()
        self.vc_plot = VCPlot()
        self.experiment_browser = ExperimentBrowser()
            # get VC traces - split into -55, -70, plot
            # get IC traces - split into -55, -70, plot
            # get IC latency from CS, plot as infinite line on all plots, fill in latency value
            # get synapse, gap from yml, fill in those parameters
            # connect fit button to fitting algorithm (which one? from where?)
                # need to pull set latency value as starting point
            # plot fit on top of traces
            # connect save button to dump result into new DB
            # connect load next pair button to clear results and move to the next pair

    def load_pair(self, pair):
        print ('loading responses...')
        s = db.Session()
        q = response_query(s)
        q = q.filter(db.PulseResponse.pair_id==pair.id)
        self.pulse_responses = q.all()
        print('got this many responses: %d' % len(self.pulse_responses))
        s.close()


    def analyze_responses(self):
        ex_limits = [-80e-3, -63e-3]
        in_limits = [-62e-3, -45e-3]
        # for rec in self.pulse_responses:

    def response_query(self, session):
        q = session.query(
        db.PulseResponse.id.label('response_id'),
        db.PulseResponse.data,
        db.PulseResponse.start_time.label('rec_start'),
        db.StimPulse.onset_time.label('pulse_start'),
        db.StimPulse.duration.label('pulse_dur'),
        db.StimSpike.max_dvdt_time.label('spike_time'),
        db.PatchClampRecording.clamp_mode,
        db.PatchClampRecording.baseline_potential,
        )
        q = q.join(db.StimPulse, db.PulseResponse.stim_pulse)
        q = q.join(db.StimSpike, db.StimSpike.stim_pulse_id==db.StimPulse.id)
        q = q.join(db.Recording, db.PulseResponse.recording)
        q = q.join(db.PatchClampRecording)




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