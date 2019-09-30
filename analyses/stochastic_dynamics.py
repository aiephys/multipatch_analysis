# coding: utf8
from __future__ import print_function, division
import sys
from collections import OrderedDict
import numpy as np
import numba
import pyqtgraph as pg
import pyqtgraph.multiprocess
from pyqtgraph.Qt import QtGui, QtCore
import scipy.stats as stats
from aisynphys.database import database as db
from aisynphys.pulse_response_strength import PulseResponseStrength
from aisynphys.synapse_prediction import SynapsePrediction, get_amps, get_baseline_amps
from aisynphys.ui.ndslicer import NDSlicer
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
        available_vesicles = int(np.clip(np.round(state['available_vesicle']), 0, self.n_release_sites))
        return release_likelihood(amplitudes, available_vesicles, self.release_probability, self.mini_amplitude, self.mini_amplitude_stdev, self.measurement_stdev)


# @numba.jit
def release_likelihood(amplitudes, available_vesicles, release_probability, mini_amplitude, mini_amplitude_stdev, measurement_stdev):
    """Return a measure of the likelihood that a synaptic response will have certain amplitude(s),
    given the state parameters for the synapse.
    
    Parameters
    ----------
    amplitudes : array
        The amplitudes for which likelihood values will be returned
    available_vesicles : int
        Number of vesicles available for release
    release_probability : float
        Probability for each available vesicle to be released
    mini_amplitude : float
        Mean amplitude of response evoked by a single vesicle release
    mini_amplitude_stdev : float
        Standard deviation of response amplitudes evoked by a single vesicle release
    measurement_stdev : float
        Standard deviation of response amplitude measurement errors
        
        
    For each value in *amplitudes*, we calculate the likelihood that a synapse would evoke a response
    of that amplitude. Likelihood is calculated as follows:
    
    1. Given the number of vesicles available to be released (nV) and the release probability (pR), determine
       the probability that each possible number of vesicles (nR) will be released using the binomial distribution
       probability mass function. For example, if there are 3 vesicles available and the release probability is
       0.1, then the possibilities are:
           vesicles released (nR)    probability
                                0    0.729
                                1    0.243
                                2    0.27
                                3    0.001
    2. For each possible number of released vesicles, calculate the likelihood that this possibility could
       evoke a response of the tested amplitude. This is calculated using the Gaussian probability distribution 
       function where µ = nR * mini_amplitude and σ = sqrt(mini_amplitude_stdev^2 * nR + measurement_stdev)
    3. The total likelihood is the sum of likelihoods for all possible values of nR.
    """
    amplitudes = np.array(amplitudes)
    
    likelihood = np.zeros(len(amplitudes))
    release_prob = stats.binom(available_vesicles, release_probability)
    
    n_vesicles = np.arange(available_vesicles + 1)
    
    # probability of releasing n_vesicles given available_vesicles and release_probability
    p_n = release_prob.pmf(n_vesicles)
    
    # expected amplitude for n_vesicles
    amp_mean = n_vesicles * mini_amplitude
    
    # amplitude stdev increases by sqrt(n) with number of released vesicles
    amp_stdev = (mini_amplitude_stdev**2 * n_vesicles + measurement_stdev**2) ** 0.5
    
    # distributions of amplitudes expected for n_vesicles
    amp_prob = p_n[None, :] * normal_pdf(amp_mean[None, :], amp_stdev[None, :], amplitudes[:, None])
    
    # sum all distributions across n_vesicles
    likelihood = amp_prob.sum(axis=1)
    
    return likelihood


def normal_pdf(mu, sigma, x):
    """Probability density function of normal distribution
    """
    return (1.0 / (2 * np.pi * sigma**2))**0.5 * np.exp(- (x-mu)**2 / (2 * sigma**2))


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
        
        self.plt4 = self.glw.addPlot(0, 1, title="amplitude distributions", rowspan=3)
        # self.plt4.setYLink(self.plt2)
        self.plt4.setMaximumWidth(500)
        self.plt4.selected_items = []

        self.amp_sample_values = np.linspace(-0.02, 0.1, 200)

    def set_result(self, model, result, pre_state, post_state):
        self.model = model
        self.result = result
        self.pre_state = pre_state
        self.post_state = post_state

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
        self.plt3.plot(compressed_spike_times, pre_state['available_vesicle'], pen=None, symbol='t', symbolBrush=brushes)
        self.plt3.plot(compressed_spike_times, post_state['available_vesicle'], pen=None, symbol='o', symbolBrush=brushes)

        self.plt4.clear()
        
        # plot full distribution of event amplitudes
        bins = np.linspace(self.amp_sample_values[0], self.amp_sample_values[-1], 40)
        amp_hist = np.histogram(self.result['amplitude'], bins=bins)
        self.plt4.plot(amp_hist[1], amp_hist[0] * len(amp_hist[0]) / amp_hist[0].sum(), stepMode=True, fillLevel=0, brush=0.3)

        # plot average model event distribution
        amps = self.amp_sample_values
        total_dist = np.zeros(len(amps))
        for i in range(self.result.shape[0]):
            state = self.pre_state[i]
            total_dist += self.model.likelihood(amps, state)
        total_dist *= len(total_dist) / total_dist.sum()
        self.plt4.plot(amps, total_dist, fillLevel=0, brush=(255, 0, 0, 50))
    
    def amp_sp_clicked(self, sp, pts):
        i = pts[0].index()
        state = self.pre_state[i]
        expected_amp = self.result[i]['expected_amplitude']
        measured_amps = self.result[i]['amplitude']
        amps = self.amp_sample_values
        
        for item in self.plt4.selected_items:
            self.plt4.removeItem(item)
        l = self.model.likelihood(amps, state)
        p = self.plt4.plot(amps, l * len(l) / l.sum(), pen=(255, 255, 0, 100))
        l1 = self.plt4.addLine(x=self.result[i]['amplitude'])
        l2 = self.plt4.addLine(x=expected_amp, pen='r')
        self.plt4.selected_items = [p, l1, l2]


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

        with pg.multiprocess.Parallelize(enumerate(all_inds), results=self.result, progressDialog='synapticulating...', workers=8) as tasker:
            for i, inds in tasker:
                params = self[inds]
                tasker.results[inds] = func(params)
        
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
            result_img[ind] = np.log(param_space.result[ind][1]['likelihood'] + 1).mean()
        self.slicer.set_data(result_img)
        self.results = result_img
        
        best = np.unravel_index(np.argmax(result_img), result_img.shape)
        self.select_result(best)
        
    def selection_changed(self, slicer):
        index = tuple(slicer.index().values())
        self.select_result(index, update_slicer=False)

    def select_result(self, index, update_slicer=True):
        result = self.param_space.result[index]
        self.result_widget.set_result(*result)
        if update_slicer:
            self.slicer.set_index(index)


if __name__ == '__main__':
    pg.mkQApp()
    pg.dbg()
    
    # strong ex, no failures, no depression
    expt_id = 1535402792.695
    pre_cell_id = 8
    post_cell_id = 7
    
    # strong ex with failures
    expt_id = 1537820585.767
    pre_cell_id = 1
    post_cell_id = 2
    
    # strong ex, depressing
    expt_id = 1536781898.381
    pre_cell_id = 8
    post_cell_id = 2

    # expt_id = float(sys.argv[1])
    # pre_cell_id = int(sys.argv[2])
    # post_cell_id = int(sys.argv[3])


    session = db.session()


    expt = db.experiment_from_timestamp(expt_id, session=session)
    pair = expt.pairs[(pre_cell_id, post_cell_id)]

    syn_type = pair.synapse.synapse_type


    # 1. Get a list of all presynaptic spike times and the amplitudes of postsynaptic responses

    raw_events = get_amps(session, pair, clamp_mode='ic')
    mask = event_qc(raw_events)
    events = raw_events[mask]

    rec_times = 1e-9 * (events['rec_start_time'].astype(float) - float(events['rec_start_time'][0]))
    spike_times = events['max_dvdt_time'] + rec_times
    amplitude_field = 'pos_dec_amp'    

    # 2. Initialize model parameters:
    #    - release model with depression, facilitation
    #    - number of synapses, distribution of per-vesicle amplitudes estimated from first pulse CV
    #    - measured distribution of background noise
    #    - parameter space to be searched

    raw_bg_events = get_baseline_amps(session, pair, clamp_mode='ic')
    mask = event_qc(raw_bg_events)
    bg_events = raw_bg_events[mask]
    mean_bg_amp = bg_events[amplitude_field].mean()
    amplitudes = events[amplitude_field] - mean_bg_amp

    first_pulse_mask = events['pulse_number'] == 1
    first_pulse_amps = amplitudes[first_pulse_mask]
    first_pulse_stdev = first_pulse_amps.std()
    

    def log_space(start, stop, steps):
        return start * (stop/start)**(np.arange(steps) / (steps-1))
        
    n_release_sites = 20
    release_probability = 0.1
    max_events = -1
    mini_amp_estimate = first_pulse_amps.mean() / (n_release_sites * release_probability)
    params = {
        'n_release_sites': np.array([1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]),
        'release_probability': log_space(0.01, 0.9, 21),
        'mini_amplitude': log_space(0.001, 0.3, 21),
        'mini_amplitude_stdev': log_space(0.0001, 0.1, 21),
        'measurement_stdev': bg_events[amplitude_field].std(),
        'recovery_tau': log_space(2e-5, 20, 21),
    }

    # quick test
    # n_release_sites = 8
    # release_probability = 0.1
    # mini_amp_estimate = first_pulse_amps.mean() / (n_release_sites * release_probability)
    # max_events = 20
    # params = {
    #     'n_release_sites': np.array([1, 2, 4, 8, 16]),
    #     'release_probability': np.array([0.1, 0.2, 0.4, 0.6, 0.8, 1.0]),
    #     'mini_amplitude': mini_amp_estimate * 1.2**np.arange(-24, 24, 2),
    #     'mini_amplitude_stdev': mini_amp_estimate * 0.2 * 1.2**np.arange(-24, 24, 8),
    #     'measurement_stdev': 0.001,
    #     'recovery_tau': 0.01,
    # }

    # # Effects of mini_amp_stdev
    # n_release_sites = 20
    # release_probability = 0.1
    # max_events = -1
    # mini_amp_estimate = first_pulse_amps.mean() / (n_release_sites * release_probability)
    # params = {
    #     'n_release_sites': np.arange(1, 30),
    #     'release_probability': release_probability,
    #     'mini_amplitude': mini_amp_estimate * 1.2**np.arange(-24, 24),
    #     'mini_amplitude_stdev': mini_amp_estimate * 1.2**np.arange(-24, 24),
    #     'measurement_stdev': 0,
    #     'recovery_tau': 0.01,
    # }

    # 3. For each point in the parameter space, simulate a synapse and estimate the joint probability of the set of measured amplitudes
    
    def run_model(params):    
        model = StochasticReleaseModel()
        for k,v in params.items():
            setattr(model, k, v)
        return (model,) + model.measure_likelihood(spike_times[:max_events], amplitudes[:max_events])
    
    
    param_space = ParameterSpace(params)
    param_space.run(run_model)

    # 4. Visualize / characterize mapped parameter space. Somehow.
    win = ParameterSearchWidget(param_space)
    win.show()
