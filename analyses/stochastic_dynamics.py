# coding: utf8
from __future__ import print_function, division
import os, sys, time, pickle, functools
from collections import OrderedDict
import numpy as np
import numba
import pyqtgraph as pg
import pyqtgraph.multiprocess
from pyqtgraph.Qt import QtGui, QtCore
import scipy.stats as stats
import scipy.optimize
from aisynphys.database import default_db as db
from aisynphys.ui.ndslicer import NDSlicer
from aisynphys import config


class StochasticReleaseModel(object):
    """
    
    Parameters
    ----------
    params : dict
        A dictionary of parameters specifying the behavior of the model:

        - n_release_sites (int) : Number of synaptic release zones
        - base_release_probability (float) : Resting-state synaptic release probability (0.0-1.0)
        - mini_amplitude (float) : Mean PSP amplitude evoked by a single vesicle release
        - mini_amplitude_stdev (float) : Stdev of PSP amplitude evoked from single releases
        - vesicle_recovery_tau (float) : Time constant for vesicle replenishment ("typically in the order of seconds" accoring to Hennig 2013)
        - facilitation_amount (float) : Release probability facilitation per spike (0.0-1.0)
        - facilitation_recovery_tau (float) : Time constant for facilitated release probability to recover toward resting state
        - measurement_stdev (float) : Extra variance in PSP amplitudes purely as a result of membrane noise / measurement error
    """
    
    result_dtype = [
        ('spike_time', float),
        ('amplitude', float),
        ('expected_amplitude', float),
        ('likelihood', float),
    ]
    
    state_dtype = [
        ('available_vesicle', float),
        ('release_probability', float),
    ]

    param_names = [
        'n_release_sites',
        'base_release_probability',
        'mini_amplitude',
        'mini_amplitude_stdev',
        'vesicle_recovery_tau',
        'facilitation_amount',
        'facilitation_recovery_tau',
        'measurement_stdev',
    ]
        
    def __init__(self, params):
        for k in params:
            if k not in self.param_names:
                raise ValueError("Unknown parameter name %r" % k)
        self.params = params

        # How long to wait after a NaN event before the model begins accumulating likelihood values again
        self.missing_event_penalty = 0.0

    def optimize_mini_amplitude(self, spike_times, amplitudes):
        """Given a set of spike times and amplitudes, optimize the mini_amplitude parameter
        to produce the highest likelihood model.
        
        Returns the output of measure_likelihood for the best model.
        """
        params = self.params.copy()

        init_amp = estimate_mini_amplitude(amplitudes, params)
        params['mini_amplitude'] = init_amp
        init_result = self.measure_likelihood(spike_times, amplitudes, params)
        mean_amp = np.nanmean(amplitudes)
        ratio = mean_amp / np.nanmean(init_result['result']['expected_amplitude'])
        init_amp *= ratio
        
        # in some cases, init_amp is way too large; force these back down:
        init_amp = min(init_amp, mean_amp) if mean_amp > 0 else max(init_amp, mean_amp)
             
        # self.params['mini_amplitude'] = params['mini_amplitude']
        result = self.optimize(spike_times, amplitudes, optimize={'mini_amplitude': (init_amp, init_amp*0.01, init_amp*100)})
        
        result['optimization_info'] = {'init_amp': init_amp / ratio, 'ratio': ratio, 'corrected_amp': init_amp, 'init_likelihood': init_result['likelihood']}
        return result
    
    def optimize(self, spike_times, amplitudes, optimize):
        """Optimize specific parameters to maximize the model likelihood.

        This method updates the attributes for any optimized parameters and returns the
        best result from measure_likelihood().
        
        Parameters
        ----------
        spike_times : array
            Times (in seconds) of presynaptic spikes in ascending order
        amplitudes : array
            Evoked PSP/PSC amplitudes for each spike listed in *spike_times*. Amplitudes may be
            NaN to indicate that the event should be ignored (usually because the spike could not be
            detected or the amplitude could not be measured). Any events within 10 seconds following
            a NaN will update the model as usual, but will not be included in the likelihood estimate.
        optimize : dict | callable
            A dictionary of {'parameter_name': (init, min, max)} that specifies model
            parameters to be optimized. Alternatively, this may be a function that accepts 
            a dictionary of fixed model parameters and returns the dictionary of parameters
            to be optimized, as described above.
        """
        params = self.params.copy()
        if callable(optimize):
            optimize = optimize(self, params)

        init = [val[0] for val in optimize.values()]
        bounds = [sorted(val[1:3]) for val in optimize.values()]
        results = {}
        
        def fn(x):
            opt_params = params.copy()
            for i,k in enumerate(optimize.keys()):
                opt_params[k] = np.clip(x[i], *bounds[i])
            res = self.measure_likelihood(spike_times, amplitudes, opt_params)
            results[tuple(x)] = res
            # print(opt_params)
            # print(res['likelihood'])
            return -res['likelihood']
        
        best = scipy.optimize.minimize(fn, x0=init, 
            #method='BFGS', bounds=bounds, options={'gtol': 1, 'eps': 10e-6}  # oscillates
            method="Nelder-Mead", options={'fatol': 0.01}  # no bounds, can't set initial step?
            #method='Powell', options={'ftol': 0.01}  # not efficient; can't set initial step
            #method='CG', options={'eps': 10e-6}  # not efficient
            #method='Newton-CG',  # requires jacobian
            #method='L-BFGS-B'  # not efficient
            #method='COBYLA', options={'rhobeg': -100e-6}  # fast but oscillates, misses peak
            
        )
        
        best_result = results[tuple(best.x.flatten())]
        best_result['optimization_init'] = optimize
        best_result['optimization_result'] = best

        # plot optimization route (for debugging)
        # x = [k[0] for k in results.keys()]
        # y = [v['likelihood'] for v in results.values()]
        # brushes = [pg.mkBrush((i, int(len(x)*1.2))) for i in range(len(x))]
        # plt = pg.plot(x, y, pen=None, symbol='o', symbolBrush=brushes)
        # plt.addLine(x=best.x[0])
        # plt.addLine(y=best_result['likelihood'])
        
        # update attributes with best result
        for i,k in enumerate(optimize.keys()):
            self.params[k] = best.x[i]
                    
        return best_result
    
    def measure_likelihood(self, spike_times, amplitudes, params=None):
        """Compute a measure of the likelihood that *times* and *amplitudes* could be generated
        by a synapse with the current dynamic parameters.
        
        Parameters
        ----------
        spike_times : array
            Times (in seconds) of presynaptic spikes in ascending order
        amplitudes : array
            Evoked PSP/PSC amplitudes for each spike listed in *spike_times*. Amplitudes may be
            NaN to indicate that the event should be ignored (usually because the spike could not be
            detected or the amplitude could not be measured). Any events within 10 seconds following
            a NaN will update the model as usual, but will not be included in the likelihood estimate.
        
        Returns
        -------
        result : array
            result contains fields: spike_time, amplitude, expected_amplitude, likelihood
        pre_spike_state:
            state variables immediately before each spike
        post_spike_state:
            state variables immediately after each spike
        params : dict
            A dictionary of parameters used to generate the model result
        """
        if params is None:
            params = self.params

        assert params['n_release_sites'] < 67, "For n_release_sites > 66 we need to use scipy.special.binom instead of the optimized binom_coeff"

        result = np.empty(len(spike_times), dtype=self.result_dtype)
        pre_spike_state = np.full(len(spike_times), np.nan, dtype=self.state_dtype)
        post_spike_state = np.full(len(spike_times), np.nan, dtype=self.state_dtype)
        
        self._run_model(
            spike_times=spike_times, 
            amplitudes=amplitudes, 
            result=result, 
            pre_spike_state=pre_spike_state, 
            post_spike_state=post_spike_state, 
            missing_event_penalty=self.missing_event_penalty,
            **params,
        )
        
        # scalar representation of overall likelihood
        likelihood = np.nanmean(np.log(result['likelihood'] + 0.1))
        
        return {
            'result': result, 
            'pre_spike_state': pre_spike_state,
            'post_spike_state': post_spike_state,
            'likelihood': likelihood,
            'params': params,
        }
            
    def likelihood(self, amplitudes, state, params=None):
        """Estimate the probability density of seeing a particular *amplitude*
        given a number of *available_vesicles*.
        """
        if params is None:
            params = self.params.copy()
        
        available_vesicles = int(np.clip(np.round(state['available_vesicle']), 0, params['n_release_sites']))
        return release_likelihood(amplitudes, available_vesicles, state['release_probability'], params['mini_amplitude'], params['mini_amplitude_stdev'], params['measurement_stdev'])


    @staticmethod
    @numba.jit(nopython=True)
    def _run_model( spike_times, 
                    amplitudes,
                    result,
                    pre_spike_state,
                    post_spike_state, 
                    missing_event_penalty,
                    n_release_sites,
                    base_release_probability,
                    mini_amplitude,
                    mini_amplitude_stdev,
                    vesicle_recovery_tau,
                    facilitation_amount,
                    facilitation_recovery_tau,
                    measurement_stdev,
                    ):

        # initialize state parameters:
        # available_vesicles is a float as a means of avoiding the need to model stochastic vesicle docking;
        # we just assume that recovery is a continuous process. (maybe we should just call this "neurotransmitter"
        # instead)
        available_vesicle = n_release_sites
        release_probability = base_release_probability
        
        previous_t = spike_times[0]
        last_nan_time = -np.inf

        for i,t in enumerate(spike_times):
            amplitude = amplitudes[i]
            if np.isnan(amplitude):
                expected_amplitude = np.nan
                likelihood = np.nan
                last_nan_time = t

            else:
                # recover vesicles up to the current timepoint
                dt = t - previous_t
                previous_t = t

                # recover vesicles up to the current timepoint
                v_recovery = np.exp(-dt / vesicle_recovery_tau)
                available_vesicle += (n_release_sites - available_vesicle) * (1.0 - v_recovery)

                # apply recovery from facilitation toward baseline release probability
                f_recovery = np.exp(-dt / facilitation_recovery_tau)
                release_probability += (base_release_probability - release_probability) * (1.0 - f_recovery)
                # prof('recover facilitation')
                
                # predict most likely amplitude for this spike (just for show)
                expected_amplitude = release_expectation_value(
                    max(0, available_vesicle),
                    release_probability,
                    mini_amplitude,
                )

                pre_available_vesicle = available_vesicle
                pre_release_probability = release_probability

                # measure likelihood of seeing this response amplitude
                av = max(0, min(n_release_sites, int(np.round(available_vesicle))))
                likelihood = release_likelihood_scalar(amplitude, av, release_probability, mini_amplitude, mini_amplitude_stdev, measurement_stdev)
                # prof('likelihood')
                assert likelihood > 0
                
                # release vesicles
                # note: we allow available_vesicle to become negative because this helps to ensure
                # that the overall likelihood will be low for such models
                depleted_vesicle = amplitude / mini_amplitude
                available_vesicle -= depleted_vesicle

                # apply spike-induced facilitation in release probability
                release_probability += (1.0 - release_probability) * facilitation_amount
                # prof('update state')
                
                assert np.isfinite(available_vesicle)
                
                # ignore likelihood for this event if it was too close to an unmeasurable response
                if t - last_nan_time < missing_event_penalty:
                    likelihood = np.nan

                # record model state immediately before spike
                pre_spike_state[i]['available_vesicle'] = pre_available_vesicle
                pre_spike_state[i]['release_probability'] = pre_release_probability
                    
                # record model state immediately after spike
                post_spike_state[i]['available_vesicle'] = available_vesicle
                post_spike_state[i]['release_probability'] = release_probability
                # prof('record')
            
            if np.isnan(available_vesicle):
                raise Exception("NaNs where they shouldn't be")
            
            # record results
            result[i]['spike_time'] = t
            result[i]['amplitude'] = amplitude
            result[i]['expected_amplitude'] = expected_amplitude
            result[i]['likelihood'] = likelihood


@numba.jit(nopython=True)
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
    return np.array([
        release_likelihood_scalar(amplitude, available_vesicles, release_probability, mini_amplitude, mini_amplitude_stdev, measurement_stdev) 
        for amplitude in amplitudes])


@numba.jit(nopython=True)
def release_likelihood_scalar(amplitude, available_vesicles, release_probability, mini_amplitude, mini_amplitude_stdev, measurement_stdev):
    """Same as release_likelihood, but optimized for a scalar amplitude argument"""
    n_vesicles = np.arange(available_vesicles + 1)
    
    # probability of releasing n_vesicles given available_vesicles and release_probability
    p_n = binom_pmf_range(available_vesicles, release_probability, available_vesicles + 1)
    
    # expected amplitude for n_vesicles
    amp_mean = n_vesicles * mini_amplitude
    
    # amplitude stdev increases by sqrt(n) with number of released vesicles
    amp_stdev = (mini_amplitude_stdev**2 * n_vesicles + measurement_stdev**2) ** 0.5
    
    # distributions of amplitudes expected for n_vesicles
    amp_prob = p_n * normal_pdf(amp_mean, amp_stdev, amplitude)
    
    # sum all distributions across n_vesicles
    likelihood = amp_prob.sum()
    
    return likelihood


# def release_distribution(available_vesicles, release_probability, mini_amplitude, mini_amplitude_stdev, measurement_stdev):
#     """Return the release amplitude distribution defined by the arguments.
#     """
#     # calculate positive and negate at the end if needed.
#     sign = mini_amplitude / abs(mini_amplitude)
#     mini_amplitude = abs(mini_amplitude)
    
#     mn = -measurement_stdev * 3
#     mx = max(0, available_vesicles) * mini_amplitude + measurement_stdev * 3
#     n_samp = int(max(0, available_vesicles) + 1) * 20
#     amplitudes = np.linspace(mn, mx, n_samp)
#     da = amplitudes[1] - amplitudes[0]
#     return amplitudes * sign, release_likelihood(amplitudes, available_vesicles, release_probability, mini_amplitude, mini_amplitude_stdev, measurement_stdev) * da

@numba.jit(nopython=True)
def release_expectation_value(available_vesicles, release_probability, mini_amplitude):
    """Return the expectation value for the release amplitude distribution defined by the arguments.
    """
    return binom_mean(available_vesicles, release_probability) * mini_amplitude
   

@numba.jit(nopython=True)
def normal_pdf(mu, sigma, x):
    """Probability density function of normal distribution
    
    Same as scipy.stats.norm(mu, sigma).pdf(x)
    """
    return (1.0 / (2 * np.pi * sigma**2))**0.5 * np.exp(- (x-mu)**2 / (2 * sigma**2))

#@functools.lru_cache(maxsize=2**14)
@numba.jit(nopython=True)
def binom_pmf_range(n, p, k):
    """Probability mass function of binomial distribution
    
    Same as scipy.stats.binom(n, p).pmf(arange(k)), but much faster.
    """
    k = np.arange(k)
    bc = np.array([binom_coeff(n,k1) for k1 in k])
    return bc * p**k * (1-p)**(n-k)


_binom_coeff_cache = np.fromfunction(scipy.special.binom, (67, 67)).astype(int)

@numba.jit(nopython=True)
def binom_coeff(n, k):
    """Binomial coefficient: n! / (k! (n-k)!)
    
    Same as scipy.special.binom, but much faster and limited to n < 67.
    """
    # note: one cold imagine writing an optimized binomial coefficient function that
    # is not limited to n < 67, but having to look out for integer overflows slows us down.
    return _binom_coeff_cache[n, k]


@numba.jit(nopython=True)
def binom_mean(n, p):
    """Expectation value of binomial distribution
    
    Same as stats.binom(n, p).mean(), but much-much faster.
    """
    return n * p


class ModelResultWidget(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)
        
        self.glw = pg.GraphicsLayoutWidget()
        self.layout.addWidget(self.glw, 0, 0)
        
        self.plt1 = self.glw.addPlot(0, 0, title="model likelihood vs compressed time")
        
        self.plt2 = self.glw.addPlot(1, 0, title="event amplitude vs compressed time")
        self.plt2.setXLink(self.plt1)
        
        self.state_key = 'release_probability'
        self.plt3 = self.glw.addPlot(2, 0, title=self.state_key + " vs compressed time")
        self.plt3.setXLink(self.plt1)
        
        self.plt4 = self.glw.addPlot(0, 1, title="amplitude distributions", rowspan=3)
        # self.plt4.setYLink(self.plt2)
        self.plt4.setMaximumWidth(500)
        self.plt4.selected_items = []

        self.amp_sample_values = np.linspace(-0.005, 0.005, 200)

    def set_result(self, result):
        self.model = result['model']
        self.result = result['result']
        self.pre_state = result['pre_spike_state']
        self.post_state = result['post_spike_state']
        self.params = result['params']

        # color events by likelihood
        cmap = pg.ColorMap([0, 1.0], [(0, 0, 0), (255, 0, 0)])
        threshold = 10
        err_colors = cmap.map((threshold - self.result['likelihood']) / threshold)
        brushes = [pg.mkBrush(c) for c in err_colors]

        # log spike intervals to make visualization a little easier
        compressed_spike_times = np.empty(len(self.result['spike_time']))
        compressed_spike_times[0] = 0.0
        np.cumsum(np.diff(self.result['spike_time'])**0.25, out=compressed_spike_times[1:])

        self.plt1.clear()
        self.plt1.plot(compressed_spike_times, self.result['likelihood'], pen=None, symbol='o', symbolBrush=brushes)
        
        self.plt2.clear()
        self.plt2.plot(compressed_spike_times, self.result['expected_amplitude'], pen=None, symbol='x', symbolPen=0.5, symbolBrush=brushes)
        amp_sp = self.plt2.plot(compressed_spike_times, self.result['amplitude'], pen=None, symbol='o', symbolBrush=brushes)
        amp_sp.scatter.sigClicked.connect(self.amp_sp_clicked)
        
        self.plt3.clear()
        self.plt3.plot(compressed_spike_times, self.pre_state[self.state_key], pen=None, symbol='t', symbolBrush=brushes)
        self.plt3.plot(compressed_spike_times, self.post_state[self.state_key], pen=None, symbol='o', symbolBrush=brushes)
        self.plt4.clear()
        
        # plot full distribution of event amplitudes
        bins = np.linspace(np.nanmin(self.result['amplitude']), np.nanmax(self.result['amplitude']), 40)
        amp_hist = np.histogram(self.result['amplitude'], bins=bins)
        self.plt4.plot(amp_hist[1], amp_hist[0] / amp_hist[0].sum(), stepMode=True, fillLevel=0, brush=0.3)

        # plot average model event distribution
        amps = self.amp_sample_values
        total_dist = np.zeros(len(amps))
        for i in range(self.result.shape[0]):
            state = self.pre_state[i]
            if not np.all(np.isfinite(tuple(state))):
                continue
            total_dist += self.model.likelihood(amps, state)
        total_dist /= total_dist.sum()
        self.plt4.plot(amps, total_dist, fillLevel=0, brush=(255, 0, 0, 50))
    
    def amp_sp_clicked(self, sp, pts):
        i = pts[0].index()
        state = self.pre_state[i]
        expected_amp = self.result[i]['expected_amplitude']
        measured_amp = self.result[i]['amplitude']
        amps = self.amp_sample_values
        
        for item in self.plt4.selected_items:
            self.plt4.removeItem(item)
        l = self.model.likelihood(amps, state)
        p = self.plt4.plot(amps, l / l.sum(), pen=(255, 255, 0, 100))
        l1 = self.plt4.addLine(x=measured_amp)
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
            result_img[ind] = param_space.result[ind]['likelihood']
        self.slicer.set_data(result_img)
        self.results = result_img
        
        best = np.unravel_index(np.argmax(result_img), result_img.shape)
        self.select_result(best)
        
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
        self.result_widget.set_result(result)
        
        print("----- Selected result: -----")
        print("  model parameters:")
        for k,v in result['params'].items():
            print("    {:30s}: {}".format(k, v))
        if 'optimization_init' in result:
            print("  initial optimization parameters:")
            for k,v in result['optimization_init'].items():
                print("    {:30s}: {}".format(k, v))
        if 'optimization_result' in result:
            opt = result['optimization_result']
            print("  optimization results:")
            print("    nfev:", opt.nfev)
            print("    message:", opt.message)
            print("    success:", opt.success)
            print("    status:", opt.status)
        if 'optimization_info' in result:
            print("  optimization info:")
            for k,v in result['optimization_info'].items():
                print("    {:30s}: {}".format(k, v))
        print("  likelihood: {}".format(result['likelihood']))
        
        if update_slicer:
            self.slicer.set_index(index)


def event_query(pair, db, session):
    q = session.query(
        db.PulseResponse,
        db.PulseResponse.ex_qc_pass,
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

    q = q.join(db.PulseResponseFit)
    q = q.join(db.StimPulse)
    q = q.join(db.Recording, db.PulseResponse.recording)
    q = q.join(db.PatchClampRecording)

    q = q.filter(db.PulseResponse.pair_id==pair.id)
    q = q.filter(db.PatchClampRecording.clamp_mode=='ic')
    
    q = q.order_by(db.Recording.start_time).order_by(db.StimPulse.onset_time)

    return q


def event_qc(events):
    mask = (events['ex_qc_pass'] == True)
    # mask = mask & (events['fit_nrmse'] < 0.6)
    return mask


def estimate_mini_amplitude(amplitudes, params):
    avg_amplitude = np.nanmean(amplitudes)
    expected = release_expectation_value(params['n_release_sites'], params['base_release_probability'], 1)
    init_amp = avg_amplitude / expected
    
    # takes care of cases where the first peak in the distribution is larger than any measured events;
    # otherwise this case is very difficult to optimize
    while abs(init_amp) > abs(avg_amplitude):
        init_amp /= 2

    return init_amp



if __name__ == '__main__':
    app = pg.mkQApp()
    if sys.flags.interactive == 1:
        pg.dbg()
    
    import argparse
    
    parser = argparse.ArgumentParser(parents=[config.parser])
    parser.add_argument('experiment_id', type=str, nargs='?')
    parser.add_argument('pre_cell_id', type=str, nargs='?')
    parser.add_argument('post_cell_id', type=str, nargs='?')
    parser.add_argument('experiment_id2', type=str, nargs='?')
    parser.add_argument('pre_cell_id2', type=str, nargs='?')
    parser.add_argument('post_cell_id2', type=str, nargs='?')
    parser.add_argument('--workers', type=int, default=None)
    parser.add_argument('--max-events', type=int, default=None, dest='max_events')
    parser.add_argument('--no-cache', default=False, action='store_true', dest='no_cache')
    
    args = parser.parse_args()
    
    
    # strong ex, no failures, no depression
    # expt_id = '1535402792.695'
    # pre_cell_id = '8'
    # post_cell_id = '7'
    
    # # strong ex with failures
    # expt_id = '1537820585.767'
    # pre_cell_id = '1'
    # post_cell_id = '2'
    
    # # strong ex, depressing
    # # expt_id = '1536781898.381'
    # # pre_cell_id = '8'
    # # post_cell_id = '2'

    # # strong in, 
    # expt_id = '1540938455.803'
    # pre_cell_id = '6' 
    # post_cell_id = '7'

    # # strong in, depressing
    # expt_id = '1530559621.966'
    # pre_cell_id = '7' 
    # post_cell_id = '6'
    
    # # ex->sst
    # expt_id = '1539987094.832'
    # pre_cell_id = '3' 
    # post_cell_id = '4'
    
    # expt_id = float(sys.argv[1])
    # pre_cell_id = int(sys.argv[2])
    # post_cell_id = int(sys.argv[3])


    expt_id = args.experiment_id
    pre_cell_id = args.pre_cell_id
    post_cell_id = args.post_cell_id

    print("Running stochastic model for %s %s %s" % (expt_id, pre_cell_id, post_cell_id))

    session = db.session()


    expt = db.experiment_from_ext_id(expt_id, session=session)
    pair = expt.pairs[(pre_cell_id, post_cell_id)]

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
    
    # 2. Initialize model parameters:
    #    - release model with depression, facilitation
    #    - number of synapses, distribution of per-vesicle amplitudes estimated from first pulse CV
    #    - measured distribution of background noise
    #    - parameter space to be searched

    amplitudes = events['dec_fit_reconv_amp'].to_numpy()
    bg_amplitudes = events['baseline_dec_fit_reconv_amp'].to_numpy()

    qc_mask = event_qc(events)
    print("%d events passed qc" % qc_mask.sum())
    amplitudes[~qc_mask] = np.nan
    amplitudes[missing_spike_mask] = np.nan
    print("%d good events to be analyzed" % np.isfinite(amplitudes).sum())

    avg_amplitude = np.nanmean(amplitudes)
    # first_pulse_mask = events['pulse_number'] == 1
    # first_pulse_amps = amplitudes[first_pulse_mask]
    # first_pulse_stdev = np.nanstd(first_pulse_amps)
    # first_pulse_mean = np.nanmean(first_pulse_amps)
    

    def log_space(start, stop, steps):
        return start * (stop/start)**(np.arange(steps) / (steps-1))


    n_release_sites = np.array([1, 2, 4, 8, 16, 32, 64])
    release_probability = np.array([0.00625, 0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0])
    search_params = {
        'n_release_sites': n_release_sites,
        'base_release_probability': release_probability,
        #'mini_amplitude': avg_amplitude * 1.2**np.arange(-12, 24, 2),  # optimized by model
        'mini_amplitude_stdev': abs(avg_amplitude) * np.array([0.05, 0.1, 0.5]),
        'measurement_stdev': np.nanstd(bg_amplitudes),
        'vesicle_recovery_tau': np.array([0.0025, 0.01, 0.04, 0.16, 0.64, 2.56]),
        'facilitation_amount': np.array([0.025, 0.05, 0.1, 0.2, 0.4]),
        'facilitation_recovery_tau': np.array([0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64]),
    }
       
    
    for k,v in search_params.items():
        if np.isscalar(v):
            assert not np.isnan(v), k
        else:
            assert not np.any(np.isnan(v)), k

    # 3. For each point in the parameter space, simulate a synapse and estimate the joint probability of the set of measured amplitudes

    trunc_spike_times = spike_times[:args.max_events]
    trunc_amplitudes = amplitudes[:args.max_events]
    
    def run_model(params):
        model = StochasticReleaseModel(params)
        result = model.optimize_mini_amplitude(trunc_spike_times, trunc_amplitudes)
        result['model'] = model
        return result


    print("Parameter space:")
    for k, v in search_params.items():
        print("   ", k, v)

    cache_file = "stochastic_cache/%s_%s_%s.pkl" % (expt_id, pre_cell_id, post_cell_id)
    if not args.no_cache and os.path.exists(cache_file):
        param_space = pickle.load(open(cache_file, 'rb'))
    else:
        print("cache miss:", cache_file)
        param_space = ParameterSpace(search_params)

        # run once to jit-precompile before measuring preformance
        run_model(param_space[(0,)*len(search_params)])

        start = time.time()
        import cProfile
        # prof = cProfile.Profile()
        # prof.enable()
        
        
        param_space.run(run_model, workers=args.workers)
        # prof.disable()
        print("Run time:", time.time() - start)
        # prof.print_stats(sort='cumulative')
   
    
        tmp = cache_file + '.tmp'
        pickle.dump(param_space, open(tmp, 'wb'))
        os.rename(tmp, cache_file)

    # if a second experiment was specified, build a combined parameter space
    if args.experiment_id2 is not None:
        param_space1 = param_space
        cache_file2 = "stochastic_cache/%s_%s_%s.pkl" % (args.experiment_id2, args.pre_cell_id2, args.post_cell_id2)
        param_space2 = pickle.load(open(cache_file2, 'rb'))
        
        params = OrderedDict()
        # params['synapse'] = [
        #     '%s_%s_%s' % (args.experiment_id, args.pre_cell_id, args.post_cell_id),
        #     '%s_%s_%s' % (args.experiment_id2, args.pre_cell_id2, args.post_cell_id2),
        # ]
        params['synapse'] = np.array([0, 1])
        params.update(param_space2.params)
        
        param_space = ParameterSpace(params)
        
        param_space.result = np.stack([param_space1.result, param_space2.result])
    
    # 4. Visualize / characterize mapped parameter space. Somehow.
        
    win = ParameterSearchWidget(param_space)
            
    # set up a few default 2D slicer views
    v1 = win.slicer.params.child('2D views').addNew()
    v1['axis 0'] = 'n_release_sites'
    v1['axis 1'] = 'base_release_probability'
    v2 = win.slicer.params.child('2D views').addNew()
    v2['axis 0'] = 'vesicle_recovery_tau'
    v2['axis 1'] = 'facilitation_recovery_tau'
    win.slicer.dockarea.moveDock(v2.viewer.dock, 'bottom', v1.viewer.dock)
    v3 = win.slicer.params.child('2D views').addNew()
    v3['axis 0'] = 'vesicle_recovery_tau'
    v3['axis 1'] = 'base_release_probability'
    v4 = win.slicer.params.child('2D views').addNew()
    v4['axis 0'] = 'facilitation_amount'
    v4['axis 1'] = 'facilitation_recovery_tau'
    win.slicer.dockarea.moveDock(v4.viewer.dock, 'bottom', v3.viewer.dock)
    
    # tur on max projection for all parameters by default
    for ch in win.slicer.params.child('max project'):
        if ch.name() == 'synapse':
            continue
        ch.setValue(True)
    
    max_like = win.results.max()
    win.slicer.histlut.setLevels(max_like * 0.95, max_like)
        
    win.show()

    if sys.flags.interactive == 0:
        app.exec_()
