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
from aisynphys.database import default_db as db
from aisynphys.ui.ndslicer import NDSlicer
from aisynphys import config


class StochasticReleaseModel(object):
    
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
    
    def __init__(self):
        # model parameters
        self.n_release_sites = 20
        self.base_release_probability = 0.5
        
        # mean and stdev of PSP amplitudes expected from a single vesicle release
        self.mini_amplitude = 50e-6
        self.mini_amplitude_stdev = 20e-6
        
        # time constant for vesicle replenishment
        # "typically in the order of seconds" accoring to Hennig 2013
        self.vesicle_recovery_tau = 500e-3
        
        # facilitation and recovery time constant
        self.facilitation_amount = 0.1
        self.facilitation_recovery_tau = 20e-3
        
        # extra variance in PSP amplitudes purely as a result of membrane noise / measurement error
        self.measurement_stdev = 100e-6
    
    def measure_likelihood(self, spike_times, amplitudes):
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
        """
        # How long to wait after a NaN event before the model begins accumulating likelihood values again
        self.missing_event_penalty = 0.0

        result = np.empty(len(spike_times), dtype=self.result_dtype)
        pre_spike_state = np.full(len(spike_times), np.nan, dtype=self.state_dtype)
        post_spike_state = np.full(len(spike_times), np.nan, dtype=self.state_dtype)
        
        
        # state parameters:
        # available_vesicles is a float as a means of avoiding the need to model stochastic vesicle docking;
        # we just assume that recovery is a continuous process. (maybe we should just call this "neurotransmitter"
        # instead)
        state = {
            'available_vesicle': self.n_release_sites,
            'release_probability': self.base_release_probability,
        }
        
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
                v_recovery = np.exp(-dt / self.vesicle_recovery_tau)
                state['available_vesicle'] += (self.n_release_sites - state['available_vesicle']) * (1.0 - v_recovery)

                # apply recovery from facilitation toward baseline release probability
                rp = state['release_probability']
                rp0 = self.base_release_probability
                f_recovery = np.exp(-dt / self.facilitation_recovery_tau)
                state['release_probability'] += (rp0 - rp) * (1.0 - f_recovery)
                
                # predict most likely amplitude for this spike (just for show)
                expected_amplitude = state['available_vesicle'] * state['release_probability'] * self.mini_amplitude
                
                # record model state immediately before spike
                for k in state:
                    pre_spike_state[i][k] = state[k]
                    
                # measure likelihood of seeing this response amplitude
                if t - last_nan_time < self.missing_event_penalty:
                    # ignore likelihood for this event if it was too close to an unmeasurable response
                    likelihood = np.nan
                else:
                    likelihood = self.likelihood([amplitude], state)[0]
                
                # release vesicles
                # note: we allow available_vesicle to become negative because this helps to ensure
                # that the overall likelihood will be low for such models
                depleted_vesicle = amplitude / self.mini_amplitude
                state['available_vesicle'] -= depleted_vesicle

                # apply spike-induced facilitation in release probability
                state['release_probability'] += (1.0 - state['release_probability']) * self.facilitation_amount

                # record model state immediately after spike
                for k in state:
                    post_spike_state[i][k] = state[k]
            
            if np.isnan(state['available_vesicle']):
                raise Exception("NaNs where they shouldn't be")
            
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
        available_vesicles = np.clip(np.round(state['available_vesicle']), 0, self.n_release_sites)
        return release_likelihood(amplitudes, int(available_vesicles), state['release_probability'], self.mini_amplitude, self.mini_amplitude_stdev, self.measurement_stdev)


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
        self.plt3.plot(compressed_spike_times, pre_state[self.state_key], pen=None, symbol='t', symbolBrush=brushes)
        self.plt3.plot(compressed_spike_times, post_state[self.state_key], pen=None, symbol='o', symbolBrush=brushes)
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
        
    def run(self, func, workers=None):
        all_inds = list(np.ndindex(self.result.shape))

        with pg.multiprocess.Parallelize(enumerate(all_inds), results=self.result, progressDialog='synapticulating...', workers=workers) as tasker:
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
            result_img[ind] = np.nanmean(np.log(param_space.result[ind][1]['likelihood'] + 1))
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


if __name__ == '__main__':
    app = pg.mkQApp()
    if sys.flags.interactive == 1:
        pg.dbg()
    
    import argparse
    
    parser = argparse.ArgumentParser(parents=[config.parser])
    parser.add_argument('experiment_id', type=str, nargs='?')
    parser.add_argument('pre_cell_id', type=str, nargs='?')
    parser.add_argument('post_cell_id', type=str, nargs='?')
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--max-events", type=int, default=None, dest='max_events')
    
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

    first_pulse_mask = events['pulse_number'] == 1
    first_pulse_amps = amplitudes[first_pulse_mask]
    first_pulse_stdev = np.nanstd(first_pulse_amps)
    

    def log_space(start, stop, steps):
        return start * (stop/start)**(np.arange(steps) / (steps-1))
        
    # n_release_sites = 20
    # release_probability = 0.1
    # max_events = -1
    # mini_amp_estimate = pair.synapse.psp_amplitude / (n_release_sites * release_probability)
    # params = {
    #     'n_release_sites': np.array([1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]),
    #     'base_release_probability': log_space(0.01, 0.9, 21),
    #     'mini_amplitude': log_space(0.001, 0.3, 21),
    #     'mini_amplitude_stdev': log_space(0.0001, 0.1, 21),
    #     'measurement_stdev': np.nanstd(bg_amplitudes),
    #     'vesicle_recovery_tau': log_space(2e-5, 20, 21),
    # }

    # quick test
    n_release_sites = 8
    release_probability = 0.1
    mini_amp_estimate = np.nanmean(first_pulse_amps) / (n_release_sites * release_probability)
    params = {
        'n_release_sites': np.array([1, 2, 4, 8, 16, 32]),
        'base_release_probability': np.array([0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]),
        'mini_amplitude': mini_amp_estimate * 1.2**np.arange(-12, 24, 2),
        'mini_amplitude_stdev': mini_amp_estimate * 0.2 * 1.2**np.arange(-12, 36, 8),
        'measurement_stdev': np.nanstd(bg_amplitudes),
        'vesicle_recovery_tau': 0.01,
    }

    for k,v in params.items():
        if np.isscalar(v):
            assert not np.isnan(v), k
        else:
            assert not np.any(np.isnan(v)), k

    # # Effects of mini_amp_stdev
    # n_release_sites = 20
    # release_probability = 0.1
    # max_events = -1
    # mini_amp_estimate = first_pulse_amps.mean() / (n_release_sites * release_probability)
    # params = {
    #     'n_release_sites': np.arange(1, 30),
    #     'base_release_probability': release_probability,
    #     'mini_amplitude': mini_amp_estimate * 1.2**np.arange(-24, 24),
    #     'mini_amplitude_stdev': mini_amp_estimate * 1.2**np.arange(-24, 24),
    #     'measurement_stdev': np.nanstd(bg_amplitudes),
    #     'vesicle_recovery_tau': 0.01,
    # }

    # 3. For each point in the parameter space, simulate a synapse and estimate the joint probability of the set of measured amplitudes
    
    def run_model(params):    
        model = StochasticReleaseModel()
        for k,v in params.items():
            if not hasattr(model, k):
                raise Exception("Invalid model parameter '%s'" % k)
            setattr(model, k, v)
        return (model,) + model.measure_likelihood(spike_times[:args.max_events], amplitudes[:args.max_events])


    print("Parameter space:")
    for k, v in params.items():
        print("   ", k, v)
    
    param_space = ParameterSpace(params)
    param_space.run(run_model, workers=args.workers)

    # 4. Visualize / characterize mapped parameter space. Somehow.
    win = ParameterSearchWidget(param_space)
    win.show()

    if sys.flags.interactive == 0:
        app.exec_()