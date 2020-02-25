# coding: utf8
from __future__ import print_function, division
import os, sys, time, pickle
from collections import OrderedDict
import numpy as np
import pyqtgraph as pg
import pyqtgraph.multiprocess
from pyqtgraph.Qt import QtGui, QtCore
# import scipy.stats as stats
# import scipy.optimize
from aisynphys.stochastic_release_model import StochasticReleaseModel
from aisynphys.database import default_db as db
from aisynphys.ui.ndslicer import NDSlicer
# from aisynphys import config


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


class ParameterSpace(object):
    def __init__(self, params):
        self.params = params
        
        static_params = {}
        for param, val in list(params.items()):
            if np.isscalar(val):
                static_params[param] = params.pop(param)
        self.static_params = static_params
        
        self.param_order = list(params.keys())
        shape = tuple([len(params[p]) for p in self.param_order])
        
        self.result = np.zeros(shape, dtype=object)
        
    def axes(self):
        return OrderedDict([(ax, {'values': self.params[ax]}) for ax in self.param_order])
        
    def run(self, func, workers=None, **kwds):
        all_inds = list(np.ndindex(self.result.shape))

        with pg.multiprocess.Parallelize(enumerate(all_inds), results=self.result, progressDialog='synapticulating...', workers=workers) as tasker:
            for i, inds in tasker:
                params = self[inds]
                tasker.results[inds] = func(params, **kwds)
        
    def __getitem__(self, inds):
        params = self.static_params.copy()
        for i,param in enumerate(self.param_order):
            params[param] = self.params[param][inds[i]]
        return params


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


def event_query(pair, db, session):
    q = session.query(
        db.PulseResponse,
        db.PulseResponse.ex_qc_pass,
        db.PulseResponse.in_qc_pass,
        db.Baseline.ex_qc_pass.label('baseline_ex_qc_pass'),
        db.Baseline.in_qc_pass.label('baseline_in_qc_pass'),
        db.PulseResponseFit.fit_amp,
        db.PulseResponseFit.dec_fit_reconv_amp,
        db.PulseResponseFit.fit_nrmse,
        db.PulseResponseFit.baseline_fit_amp,
        db.PulseResponseFit.baseline_dec_fit_reconv_amp,
        db.StimPulse.first_spike_time,
        db.StimPulse.pulse_number,
        db.StimPulse.onset_time,
        db.Recording.start_time.label('rec_start_time'),
        db.PatchClampRecording.baseline_current,
    )
    q = q.join(db.Baseline, db.PulseResponse.baseline)
    q = q.join(db.PulseResponseFit)
    q = q.join(db.StimPulse)
    q = q.join(db.Recording, db.PulseResponse.recording)
    q = q.join(db.PatchClampRecording)

    q = q.filter(db.PulseResponse.pair_id==pair.id)
    q = q.filter(db.PatchClampRecording.clamp_mode=='ic')
    
    q = q.order_by(db.Recording.start_time).order_by(db.StimPulse.onset_time)

    return q


class StochasticModelRunner:
    """Handles loading data for a synapse and executing the model across a parameter space.
    """
    def __init__(self, experiment_id, pre_cell_id, post_cell_id, workers=None):
        self.experiment_id = experiment_id
        self.pre_cell_id = pre_cell_id
        self.post_cell_id = post_cell_id
        self.title = "%s %s %s" % (experiment_id, pre_cell_id, post_cell_id)
        
        self.workers = workers
        self.max_events = None
        
        self._synapse_events = None
        self._parameters = None
        self._param_space = None

    @property
    def param_space(self):
        """A ParameterSpace instance containing the model output over the entire parameter space.
        """
        if self._param_space is None:
            self._param_space = self._generate_param_space()
        return self._param_space

    def _generate_param_space(self):
        search_params = self.parameters
        
        param_space = ParameterSpace(search_params)

        # run once to jit-precompile before measuring preformance
        self.run_model(param_space[(0,)*len(search_params)])

        start = time.time()
        import cProfile
        # prof = cProfile.Profile()
        # prof.enable()
        
        param_space.run(self.run_model, workers=self.workers)
        # prof.disable()
        print("Run time:", time.time() - start)
        # prof.print_stats(sort='cumulative')
        
        return param_space

    def store_result(self, cache_file):    
        tmp = cache_file + '.tmp'
        pickle.dump(self.param_space, open(tmp, 'wb'))
        os.rename(tmp, cache_file)

    def load_result(self, cache_file):
        self._param_space = pickle.load(open(cache_file, 'rb'))
        
    @property
    def synapse_events(self):
        """Tuple containing (spike_times, amplitudes, baseline_amps)
        """
        if self._synapse_events is None:
            self._synapse_events = self._load_synapse_events()
        return self._synapse_events

    def _load_synapse_events(self):
        session = db.session()

        expt = db.experiment_from_ext_id(self.experiment_id, session=session)
        pair = expt.pairs[(self.pre_cell_id, self.post_cell_id)]

        syn_type = pair.synapse.synapse_type
        print("Synapse type:", syn_type)

        # 1. Get a list of all presynaptic spike times and the amplitudes of postsynaptic responses

        events = event_query(pair, db, session).dataframe()
        print("loaded %d events" % len(events))

        rec_times = (events['rec_start_time'] - events['rec_start_time'].iloc[0]).dt.total_seconds().to_numpy()
        spike_times = events['first_spike_time'].to_numpy() + rec_times

        # any missing spike times get filled in with the average latency
        missing_spike_mask = np.isnan(spike_times)
        print("%d events missing spike times" % missing_spike_mask.sum())
        avg_spike_latency = np.nanmedian(events['first_spike_time'] - events['onset_time'])
        pulse_times = events['onset_time'] + avg_spike_latency + rec_times
        spike_times[missing_spike_mask] = pulse_times[missing_spike_mask]

        # get individual event amplitudes
        amplitudes = events['dec_fit_reconv_amp'].to_numpy()
        
        # filter events by inhibitory or excitatory qc
        qc_field = syn_type + '_qc_pass'
        qc_mask = events[qc_field] == True
        print("%d events passed qc" % qc_mask.sum())
        amplitudes[~qc_mask] = np.nan
        amplitudes[missing_spike_mask] = np.nan
        print("%d good events to be analyzed" % np.isfinite(amplitudes).sum())

        # get background events for determining measurement noise
        bg_amplitudes = events['baseline_dec_fit_reconv_amp'].to_numpy()
        # filter by qc
        bg_qc_mask = events['baseline_'+qc_field] == True
        bg_amplitudes[~qc_mask] = np.nan
        
        # first_pulse_mask = events['pulse_number'] == 1
        # first_pulse_amps = amplitudes[first_pulse_mask]
        # first_pulse_stdev = np.nanstd(first_pulse_amps)
        # first_pulse_mean = np.nanmean(first_pulse_amps)
        
        if self.max_events is not None:
            spike_times = spike_times[:self.max_events]
            amplitudes = amplitudes[:self.max_events]
        
        return spike_times, amplitudes, bg_amplitudes

    @property
    def parameters(self):
        """A structure defining the parameters to search.
        """
        if self._parameters is None:
            self._parameters = self._generate_parameters()
        return self._parameters  

    def _generate_parameters(self):
        spike_times, amplitudes, bg_amplitudes = self.synapse_events
        
        search_params = {
            'n_release_sites': np.array([1, 2, 4, 8, 16, 32, 64]),
            #'n_release_sites': np.array([1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]),
            'base_release_probability': np.array([0.00625, 0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]),
            #'base_release_probability': 1.0 / 1.5**(np.arange(15)[::-1]),
            #'mini_amplitude': avg_amplitude * 1.2**np.arange(-12, 24, 2),  # optimized by model
            'mini_amplitude_cv': np.array([0.05, 0.1, 0.2, 0.4, 0.8]),
            'measurement_stdev': np.nanstd(bg_amplitudes),
            'vesicle_recovery_tau': np.array([0.0025, 0.01, 0.04, 0.16, 0.64, 2.56]),
            'facilitation_amount': np.array([0.0, 0.00625, 0.025, 0.05, 0.1, 0.2, 0.4]),
            'facilitation_recovery_tau': np.array([0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28]),
        }
        
        # sanity checking
        for k,v in search_params.items():
            if np.isscalar(v):
                assert not np.isnan(v), k
            else:
                assert not np.any(np.isnan(v)), k

        print("Parameter space:")
        for k, v in search_params.items():
            print("   ", k, v)

        return search_params

    def run_model(self, params, full_result=False, **kwds):
        model = StochasticReleaseModel(params)
        spike_times, amplitudes, bg = self.synapse_events
        if 'mini_amplitude' in params:
            result = model.measure_likelihood(spike_times, amplitudes, **kwds)
        else:
            result = model.optimize_mini_amplitude(spike_times, amplitudes, **kwds)
        if full_result:
            result['model'] = model
            return result
        else:
            return {'likelihood': result['likelihood'], 'params': result['params']}


class CombinedModelRunner:
    """Model runner combining the results from multiple StochasticModelRunner instances.
    """
    def __init__(self, runners):
        self.model_runners = runners
        self.title = " : ".join(r.title for r in runners)
        
        params = OrderedDict()
        # params['synapse'] = [
        #     '%s_%s_%s' % (args.experiment_id, args.pre_cell_id, args.post_cell_id),
        #     '%s_%s_%s' % (args.experiment_id2, args.pre_cell_id2, args.post_cell_id2),
        # ]
        params['synapse'] = np.arange(len(runners))
        params.update(runners[0].param_space.params)
        param_space = ParameterSpace(params)
        param_space.result = np.stack([runner.param_space.result for runner in runners])
        
        self.param_space = param_space
        
    def run_model(self, params, full_result=False):
        params = params.copy()
        runner = self.model_runners[params.pop('synapse')]
        return runner.run_model(params, full_result)
