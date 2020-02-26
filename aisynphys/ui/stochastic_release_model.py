# coding: utf8
from __future__ import print_function, division
import numpy as np
import pyqtgraph as pg
import pyqtgraph.multiprocess
from pyqtgraph.Qt import QtGui, QtCore
from aisynphys.stochastic_release_model import StochasticReleaseModel
from aisynphys.ui.ndslicer import NDSlicer


class ModelDisplayWidget(QtGui.QWidget):
    """UI containing an NDSlicer for visualizing the complete model parameter space, and
    a ModelSingleResultWidget for showing more detailed results from individual points in the parameter space.
    """
    def __init__(self, model_runner):
        QtGui.QWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
        self.splitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.layout.addWidget(self.splitter)
        
        self.slicer = NDSlicer(model_runner.param_space.axes())
        self.slicer.selection_changed.connect(self.selection_changed)
        self.splitter.addWidget(self.slicer)
        
        self.result_widget = ModelSingleResultWidget()
        self.splitter.addWidget(self.result_widget)

        # set up a few default 2D slicer views
        v1 = self.slicer.params.child('2D views').addNew()
        v1['axis 0'] = 'n_release_sites'
        v1['axis 1'] = 'base_release_probability'
        v2 = self.slicer.params.child('2D views').addNew()
        v2['axis 0'] = 'vesicle_recovery_tau'
        v2['axis 1'] = 'facilitation_recovery_tau'
        self.slicer.dockarea.moveDock(v2.viewer.dock, 'bottom', v1.viewer.dock)
        v3 = self.slicer.params.child('2D views').addNew()
        v3['axis 0'] = 'vesicle_recovery_tau'
        v3['axis 1'] = 'base_release_probability'
        v4 = self.slicer.params.child('2D views').addNew()
        v4['axis 0'] = 'facilitation_amount'
        v4['axis 1'] = 'facilitation_recovery_tau'
        self.slicer.dockarea.moveDock(v4.viewer.dock, 'bottom', v3.viewer.dock)
        
        # turn on max projection for all parameters by default
        for ch in self.slicer.params.child('max project'):
            if ch.name() == 'synapse':
                continue
            ch.setValue(True)
        
        self.model_runner = model_runner
        self.param_space = model_runner.param_space
        
        result_img = np.zeros(self.param_space.result.shape)
        for ind in np.ndindex(result_img.shape):
            result_img[ind] = self.param_space.result[ind]['likelihood']
        self.slicer.set_data(result_img)
        self.results = result_img
        
        # select best result
        best = np.unravel_index(np.argmax(result_img), result_img.shape)
        self.select_result(best)

        max_like = self.results.max()
        
        # if results are combined across synapses, set up colors
        if 'synapse' in self.param_space.params:
            self.slicer.params['color axis', 'axis'] = 'synapse'
            syn_axis = list(self.param_space.params.keys()).index('synapse')
            max_img = np.array([result_img.take(i, axis=syn_axis).max() for i in range(result_img.shape[syn_axis])])
            max_like = max_img.min()
            max_img = max_img.min() / max_img
            syns = self.param_space.params['synapse']
            for i in syns:
                c = pg.colorTuple(pg.intColor(i, len(syns)*1.2))
                c = pg.mkColor(c[0]*max_img[i], c[1]*max_img[i], c[2]*max_img[i])
                self.slicer.params['color axis', 'colors', str(i)] = c
            
        # set histogram range
        self.slicer.histlut.setLevels(max_like * 0.85, max_like)

    def selection_changed(self, slicer):
        index = self.selected_index()
        self.select_result(index, update_slicer=False)

    def selected_index(self):
        return tuple(self.slicer.index().values())

    def get_result(self, index=None):
        if index is None:
            index = self.selected_index()
        return self.param_space.result[index]

    def select_result(self, index, update_slicer=True):
        result = self.get_result(index)
        result['params'].update(self.param_space[index])
        
        # re-run the model to get the complete results
        full_result = self.model_runner.run_model(result['params'], full_result=True, show=True)
        self.result_widget.set_result(full_result)
        
        print("----- Selected result: -----")
        print("  model parameters:")
        for k,v in full_result['params'].items():
            print("    {:30s}: {}".format(k, v))
        if 'optimized_params' in full_result:
            print("  optimized parameters:")
            for k,v in full_result['optimized_params'].items():
                print("    {:30s}: {}".format(k, v))
        if 'optimization_init' in full_result:
            print("  initial optimization parameters:")
            for k,v in full_result['optimization_init'].items():
                print("    {:30s}: {}".format(k, v))
        if 'optimization_result' in full_result:
            opt = full_result['optimization_result']
            print("  optimization results:")
            print("    nfev:", opt.nfev)
            print("    message:", opt.message)
            print("    success:", opt.success)
            print("    status:", opt.status)
        if 'optimization_info' in full_result:
            print("  optimization info:")
            for k,v in full_result['optimization_info'].items():
                print("    {:30s}: {}".format(k, v))
        print("  likelihood: {}".format(full_result['likelihood']))
        
        if update_slicer:
            self.slicer.set_index(index)


class ModelSingleResultWidget(QtGui.QWidget):
    """Plots event amplitudes and distributions for a single stochastic model run.
    """
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)

        self.splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.splitter)
        
        self.btn_layout = QtGui.QVBoxLayout()
        self.btn_widget = QtGui.QWidget()
        self.btn_widget.setLayout(self.btn_layout)
        self.splitter.addWidget(self.btn_widget)
        
        self.opt_path_plot = pg.PlotWidget()

        self.panels = {
            'events': ModelEventPlot(self),
            'correlation': ModelEventCorrelationPlot(self),
            'optimization': ModelOptimizationPlot(self),
        }
        for name, panel in self.panels.items():
            btn = QtGui.QPushButton(name)
            self.btn_layout.addWidget(btn)
            btn.setCheckable(True)
            btn.toggled.connect(panel.set_visible)
            btn.setMaximumWidth(150)
            panel.show_btn = btn
        
        self.panels['events'].show_btn.setChecked(True)

    def set_result(self, result):
        self.result = result
        for p in self.panels.values():
            p.result_changed()


class ModelResultView(object):
    """Displays one aspect of a model result.
    """
    def __init__(self, parent, visible=False):
        """Responsible for attaching widgets to parent.
        """
        self._parent = parent
        self._visible = False
        self.widget = QtGui.QWidget()
        self.layout = QtGui.QGridLayout()
        self.widget.setLayout(self.layout)
        parent.splitter.addWidget(self.widget)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self._need_update = False
        
        self.set_visible(visible)
        
    def set_visible(self, visible):
        self._visible = visible
        self.widget.setVisible(visible)
        if visible and self._need_update:
            self.update_display()

    @property
    def visible(self):
        return self._visible

    def result_changed(self):
        """Called when the model result has been set on parent
        """
        self._need_update = True
        if self.visible:
            self.update_display()
        
    def update_display(self):
        """Extend in subclass to update displayed results (taken from self._parent.result)
        """
        self._need_update = False    
        

class ModelEventPlot(ModelResultView):
    def __init__(self, parent):
        ModelResultView.__init__(self, parent)

        self.event_view = pg.GraphicsLayoutWidget()
        self.layout.addWidget(self.event_view)
        
        self.plt1 = self.event_view.addPlot(0, 0, title="model likelihood vs compressed time")
        
        self.plt2 = self.event_view.addPlot(1, 0, title="event amplitude vs compressed time")
        self.plt2.setXLink(self.plt1)
        
        self.state_key = 'release_probability'
        self.plt3 = self.event_view.addPlot(2, 0, title=self.state_key + " vs compressed time")
        self.plt3.setXLink(self.plt1)
        
        self.plt4 = self.event_view.addPlot(0, 1, title="amplitude distributions", rowspan=3)
        # self.plt4.setYLink(self.plt2)
        self.plt4.setMaximumWidth(500)
        self.plt4.selected_items = []

        self.amp_sample_values = np.linspace(-0.005, 0.005, 800)

    def update_display(self):
        ModelResultView.update_display(self)

        full_result = self._parent.result
        model = full_result['model']
        result = full_result['result']
        pre_state = full_result['pre_spike_state']
        post_state = full_result['post_spike_state']
        
        # color events by likelihood
        cmap = pg.ColorMap([0, 1.0], [(0, 0, 0), (255, 0, 0)])
        threshold = 10
        err_colors = cmap.map((threshold - result['likelihood']) / threshold)
        brushes = [pg.mkBrush(c) for c in err_colors]

        # log spike intervals to make visualization a little easier
        compressed_spike_times = np.empty(len(result['spike_time']))
        compressed_spike_times[0] = 0.0
        np.cumsum(np.diff(result['spike_time'])**0.25, out=compressed_spike_times[1:])

        self.plt1.clear()
        self.plt1.plot(compressed_spike_times, result['likelihood'], pen=None, symbol='o', symbolBrush=brushes)
        
        self.plt2.clear()
        self.plt2.plot(compressed_spike_times, result['expected_amplitude'], pen=None, symbol='x', symbolPen=0.5, symbolBrush=brushes)
        amp_sp = self.plt2.plot(compressed_spike_times, result['amplitude'], pen=None, symbol='o', symbolBrush=brushes)
        amp_sp.scatter.sigClicked.connect(self.amp_sp_clicked)
        
        self.plt3.clear()
        self.plt3.plot(compressed_spike_times, pre_state[self.state_key], pen=None, symbol='t', symbolBrush=brushes)
        self.plt3.plot(compressed_spike_times, post_state[self.state_key], pen=None, symbol='o', symbolBrush=brushes)
        self.plt4.clear()
        
        # plot full distribution of event amplitudes
        bins = np.linspace(np.nanmin(result['amplitude']), np.nanmax(result['amplitude']), 40)
        d_amp = bins[1] - bins[0]
        amp_hist = np.histogram(result['amplitude'], bins=bins)
        self.plt4.plot(amp_hist[1], amp_hist[0] / (amp_hist[0].sum() * d_amp), stepMode=True, fillLevel=0, brush=0.3)

        # plot average model event distribution
        amps = self.amp_sample_values
        d_amp = amps[1] - amps[0]
        total_dist = np.zeros(len(amps))
        for i in range(result.shape[0]):
            state = pre_state[i]
            if not np.all(np.isfinite(tuple(state))):
                continue
            total_dist += model.likelihood(amps, state)
        total_dist /= total_dist.sum() * d_amp
        self.plt4.plot(amps, total_dist, fillLevel=0, brush=(255, 0, 0, 50))
    
    def amp_sp_clicked(self, sp, pts):
        result = self._parent.result['result']
        i = pts[0].index()
        state = self.pre_state[i]
        expected_amp = result[i]['expected_amplitude']
        measured_amp = result[i]['amplitude']
        amps = self.amp_sample_values
        
        for item in self.plt4.selected_items:
            self.plt4.removeItem(item)
        l = self.model.likelihood(amps, state)
        p = self.plt4.plot(amps, l / l.sum(), pen=(255, 255, 0, 100))
        l1 = self.plt4.addLine(x=measured_amp)
        l2 = self.plt4.addLine(x=expected_amp, pen='r')
        self.plt4.selected_items = [p, l1, l2]


class ModelOptimizationPlot(ModelResultView):
    """Plots the mini amplitude optimization path
    """
    def __init__(self, parent):
        ModelResultView.__init__(self, parent)
        self.plot = pg.PlotWidget()
        self.layout.addWidget(self.plot)
        
    def update_display(self):
        ModelResultView.update_display(self)

        result = self._parent.result
        plt = self.plot
        x = result['optimization_path']['mini_amplitude']
        y = result['optimization_path']['likelihood']
        brushes = [pg.mkBrush((i, int(len(x)*1.2))) for i in range(len(x))]
        plt.clear()
        plt.plot(x, y, pen=None, symbol='o', symbolBrush=brushes)
        plt.addLine(x=result['optimized_params']['mini_amplitude'])
        plt.addLine(y=result['likelihood'])


class ModelEventCorrelationPlot(ModelResultView):
    def __init__(self, parent):
        ModelResultView.__init__(self, parent)
        self.plot = pg.PlotWidget(labels={'left': 'amp 2', 'bottom': 'amp 1'})
        self.plot.showGrid(True, True)
        self.layout.addWidget(self.plot)
        
    def update_display(self):
        ModelResultView.update_display(self)
        result = self._parent.result
        spikes = result['result']['spike_time']
        amps = result['result']['amplitude']
        
        cmap = pg.ColorMap(np.linspace(0, 1, 4), np.array([[255, 255, 255], [255, 255, 0], [255, 0, 0], [0, 0, 0]], dtype='ubyte'))
        cvals = [((np.log(dt) / np.log(10)) + 2) / 4. for dt in np.diff(spikes)]
        brushes = [pg.mkBrush(cmap.map(c)) for c in cvals]
        
        self.plot.clear()
        self.plot.plot(amps[:-1], amps[1:], pen=None, symbol='o', symbolBrush=brushes)
