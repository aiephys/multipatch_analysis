from __future__ import print_function, division
import sys
from collections import OrderedDict
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import scipy.stats as stats
from multipatch_analysis.database import database as db
from multipatch_analysis.pulse_response_strength import PulseResponseStrength
from multipatch_analysis.connection_strength import ConnectionStrength, get_amps, get_baseline_amps
from multipatch_analysis.ui.ndslicer import NDSlicer
from neuroanalysis.synaptic_release import ReleaseModel


class StochasticReleaseModel(object):
    
    result_dtype = [
        ('spike_time', float),
        ('amplitude', float),
        ('expected_amplitude', float),
        ('likelihood', float),
    ]
    
    state_dtype = [
        ('available_vesicle', float),
    ]
    
    def __init__(self):
        # model parameters
        self.n_release_sites = 20
        self.release_probability = 0.1
        self.mini_amplitude = 50e-6
        self.mini_amplitude_stdev = 20e-6
        self.recovery_tau = 100e-3
        self.measurement_stdev = 100e-6
    
    def measure_likelihood(self, spike_times, amplitudes):
        """Compute a measure of the likelihood that *times* and *amplitudes* could be generated
        by a synapse with the current dynamic parameters.
        
        Returns
        -------
        result : array
            result contains fields: spike_time, amplitude, expected_amplitude, likelihood
        pre_spike_state:
            state variables immediately before each spike
        post_spike_state:
            state variables immediately after each spike
        """
        result = np.empty(len(spike_times), dtype=self.result_dtype)
        pre_spike_state = np.empty(len(spike_times), dtype=self.state_dtype)
        post_spike_state = np.empty(len(spike_times), dtype=self.state_dtype)
        
        # state parameters:
        # available_vesicles is a float as a means of avoiding the need to model stochastic vesicle docking;
        # we just assume that recovery is a continuous process. (maybe we should just call this "neurotransmitter"
        # instead)
        state = {
            'available_vesicle': self.n_release_sites,
        }
        
        previous_t = spike_times[0]
        
        for i,t in enumerate(spike_times):
            amplitude = amplitudes[i]

            # recover vesicles up to the current timepoint
            dt = t - previous_t
            previous_t = t
            recovery = np.exp(-dt / self.recovery_tau)
            state['available_vesicle'] = state['available_vesicle'] * recovery + self.n_release_sites * (1.0 - recovery)
            
            # record model state immediately before spike
            for k in state:
                pre_spike_state[i][k] = state[k]
                
            # measure likelihood of seeing this response amplitude
            expected_amplitude = state['available_vesicle'] * self.release_probability * self.mini_amplitude
            likelihood = self.likelihood([amplitude], state)[0]
            
            # release vesicles
            # note: we allow available_vesicle to become negative because this help to ensure
            # that the overall likelihood will be low for such models
            depleted_vesicle = amplitude / self.mini_amplitude
            state['available_vesicle'] -= depleted_vesicle

            # record model state immediately after spike
            for k in state:
                post_spike_state[i][k] = state[k]
            
            # record results
            result[i]['spike_time'] = t
            result[i]['amplitude'] = amplitude
            result[i]['expected_amplitude'] = expected_amplitude
            result[i]['likelihood'] = likelihood
            
        return result, pre_spike_state, post_spike_state
            
    def likelihood(self, amplitudes, state):
        """Estimate the probability density of seeing a particular *amplitude*
        given a number of *available_vesicles*.
        """
        available_vesicles = int(max(0, state['available_vesicle']))
        
        likelihood = np.zeros(len(amplitudes))
        release_prob = stats.binom(available_vesicles, self.release_probability)
        for n_vesicles in range(available_vesicles):
            # probability of releasing n_vesicles given available_vesicles and release_probability
            p_n = release_prob.pmf(n_vesicles)
            
            # expected amplitude for n_vesicles
            amp_mean = n_vesicles * self.mini_amplitude
            
            # distribution of amplitudes expected for n_vesicles
            amp_stdev = (self.mini_amplitude_stdev**2 + self.measurement_stdev**2) ** 0.5
            amp_prob = p_n * stats.norm(amp_mean, amp_stdev).pdf(amplitudes)
            
            # increment likelihood of seeing this amplitude
            likelihood += amp_prob
        
        return likelihood


def event_qc(events):
    mask = events['ex_qc_pass'] == True
    # need more stringent qc for dynamics:
    mask &= np.abs(events['baseline_current']) < 500e-12
    return mask   


class ModelResultWidget(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
        
        self.glw = pg.GraphicsLayoutWidget()
        self.layout.addWidget(self.glw, 0, 0)
        
        self.plt1 = self.glw.addPlot(0, 0, title="likelihood vs compressed time")
        
        self.plt2 = self.glw.addPlot(1, 0, title="deconvolved amplitude vs compressed time")
        self.plt2.setXLink(self.plt1)
        
        self.plt3 = self.glw.addPlot(2, 0, title="available_vesicles vs compressed time")
        self.plt3.setXLink(self.plt1)
        
        self.plt4 = self.glw.addPlot(1, 1, title="expected amplitude distribution")
        self.plt4.setYLink(self.plt2)
        self.plt4.setMaximumWidth(300)

    def set_result(self, model, result, pre_state, post_state):
        self.model = model
        self.result = result
        self.pre_state = pre_state
        self.post_state = post_state

        # color events by likelihood
        cmap = pg.ColorMap([0, 1.0], [(0, 0, 0), (255, 0, 0)])
        err_colors = cmap.map((10 - result['likelihood']) / 10.)
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
        self.plt3.plot(compressed_spike_times, pre_state['available_vesicle'], pen=None, symbol='t', symbolBrush=brushes)
        self.plt3.plot(compressed_spike_times, post_state['available_vesicle'], pen=None, symbol='o', symbolBrush=brushes)
    
    def amp_sp_clicked(self, sp, pts):
        i = pts[0].index()
        state = self.pre_state[i]
        expected_amp = self.result[i]['expected_amplitude']
        amps = np.linspace(-expected_amp, expected_amp * 4, 2000)
        self.plt4.clear()
        self.plt4.plot(self.model.likelihood(amps, state), amps)
        self.plt4.addLine(y=self.result[i]['amplitude'])
        self.plt4.addLine(y=expected_amp, pen='r')


class PixelSelector(QtGui.QGraphicsRectItem):
    
    sigPixelSelectionChanged = QtCore.Signal(object, object)  # self, (x, y)
    
    def __init__(self, image=None, pen='y'):
        self.image = None
        QtGui.QGraphicsRectItem.__init__(self, QtCore.QRectF(0, 0, 1, 1))
        self.setImage(image)
        self.setPen(pen)
        
    def selectedPos(self):
        """Return the currently selected data location (row, col).
        """
        if self.image is None or self.image.width() == 0 or self.image.height() == 0:
            return (np.nan, np.nan)
        dataPos = self.image.mapToData(self.pos())
        return (dataPos.y(), dataPos.x())
        
    def setPen(self, pen):
        QtGui.QGraphicsRectItem.setPen(self, pg.mkPen(pen))
        
    def setImage(self, image):
        if self.scene() is not None:
            self.scene().sigMouseClicked.disconnect(self._sceneClicked)
        self.image = image
        if image is not None:
            self.setParentItem(image)
            self.scene().sigMouseClicked.connect(self._sceneClicked)            
        self.imageChanged()
        
    def imageChanged(self):
        # check new image bounds
        if self.image is None or self.image.width() == 0 or self.image.height() == 0:
            pos = (np.nan, np.nan)
        else:
            pos = [min(self.pos().x(), self.image.width()), min(self.pos.y(), self.image.height())]
        
        self.setPos(*pos)
        
    def setPos(self, *args):
        prevPos = self.pos()
        QtGui.QGraphicsRectItem.setPos(*args)
        if self.pos() != prevPos():
            self.sigPixelSelectionChanged.emit(self, self.selectedPixel())

    def _sceneClicked(self, event):
        spos = event.scenePos()
        imgPos = self.image.mapFromScene(spos)
        i, j = int(imgPos.x()), int(imgPos.y())
        self.setPos(i, j)

        
class ParameterSpace(object):
    def __init__(self, params):
        self.params = params
        
        static_params = {}
        for param, val in params.items():
            if np.isscalar(val):
                static_params[param] = params.pop(param)
        self.static_params = static_params
        
        self.param_order = list(params.keys())
        shape = tuple([len(params[p]) for p in self.param_order])
        
        self.result = np.zeros(shape, dtype=object)
        
    def axes(self):
        return OrderedDict([(ax, {'values': self.params[ax]}) for ax in self.param_order])
        
    def run(self, func):
        all_inds = list(np.ndindex(self.result.shape))
        for i,inds in enumerate(all_inds):
            params = self[inds]
            self.result[inds] = func(params)
            yield i, inds, len(all_inds)
        
    def __getitem__(self, inds):
        params = self.static_params.copy()
        for i,param in enumerate(self.param_order):
            params[param] = self.params[param][inds[i]]
        return params


class ParameterSearchWidget(QtGui.QWidget):
    def __init__(self, param_space):
        QtGui.QWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
        self.splitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.layout.addWidget(self.splitter)
        
        self.slicer = NDSlicer(param_space.axes())
        self.slicer.selection_changed.connect(self.selection_changed)
        self.splitter.addWidget(self.slicer)
        
        self.result_widget = ModelResultWidget()
        self.splitter.addWidget(self.result_widget)
        
        self.param_space = param_space
        
        result_img = np.zeros(param_space.result.shape)
        for ind in np.ndindex(result_img.shape):
            result_img[ind] = param_space.result[ind][1]['likelihood'].mean()
        self.slicer.set_data(result_img)
        
    def selection_changed(self, slicer):
        index = slicer.index()
        self.select_result(*index.values())

    def select_result(self, i, j):
        result = self.param_space.result[i, j]
        self.result_widget.set_result(*result)


if __name__ == '__main__':
    pg.mkQApp()
    pg.dbg()
    
    expt_id = 1535402792.695
    pre_cell_id = 8
    post_cell_id = 7

    # expt_id = float(sys.argv[1])
    # pre_cell_id = int(sys.argv[2])
    # post_cell_id = int(sys.argv[3])


    session = db.Session()


    expt = db.experiment_from_timestamp(expt_id, session=session)
    pair = expt.pairs[(pre_cell_id, post_cell_id)]

    syn_type = pair.connection_strength.synapse_type


    # 1. Get a list of all presynaptic spike times and the amplitudes of postsynaptic responses

    raw_events = get_amps(session, pair, clamp_mode='ic')
    mask = event_qc(raw_events)
    events = raw_events[mask]

    rec_times = 1e-9 * (events['rec_start_time'].astype(float) - float(events['rec_start_time'][0]))
    spike_times = events['max_dvdt_time'] + rec_times
    amplitude_field = 'pos_dec_amp'
    amplitudes = events[amplitude_field]
    

    # 2. Initialize model parameters:
    #    - release model with depression, facilitation
    #    - number of synapses, distribution of per-vesicle amplitudes estimated from first pulse CV
    #    - measured distribution of background noise
    #    - parameter space to be searched

    first_pulse_mask = events['pulse_number'] == 1
    first_pulse_amps = amplitudes[first_pulse_mask]
    first_pulse_stdev = first_pulse_amps.std()

    raw_bg_events = get_baseline_amps(session, pair, clamp_mode='ic')
    mask = event_qc(raw_bg_events)
    bg_events = raw_bg_events[mask]

    def log_space(start, stop, steps):
        return start * (stop/start)**(np.arange(steps) / (steps-1))
    n_release_sites = 20
    release_probability = 0.1
    mini_amp_estimate = first_pulse_amps.mean() / (n_release_sites * release_probability)
    params = {
        'n_release_sites': np.array([1, 2, 3, 4, 6, 8, 12, 16, 24, 32]),
        'release_probability': release_probability,
        'mini_amplitude': mini_amp_estimate,
        'mini_amplitude_stdev': mini_amp_estimate / 3.,
        'measurement_stdev': bg_events[amplitude_field].std(),
        'recovery_tau': log_space(2e-5, 20, 11),
    }


    # 3. For each point in the parameter space, simulate a synapse and estimate the joint probability of the set of measured amplitudes
    
    def run_model(params):    
        model = StochasticReleaseModel()
        for k,v in params.items():
            setattr(model, k, v)
        return (model,) + model.measure_likelihood(spike_times[:30], amplitudes[:30])
    
    
    param_space = ParameterSpace(params)
    for i, inds, n_inds in param_space.run(run_model):
        print(i+1, n_inds)

    # 4. Visualize / characterize mapped parameter space. Somehow.
    win = ParameterSearchWidget(param_space)
    win.show()
