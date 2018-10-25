import sys
import numpy as np
import pyqtgraph as pg
import scipy.stats as stats
from multipatch_analysis.database import database as db
from multipatch_analysis.pulse_response_strength import PulseResponseStrength
from multipatch_analysis.connection_strength import ConnectionStrength, get_amps, get_baseline_amps
from neuroanalysis.synaptic_release import ReleaseModel


class StochasticReleaseModel(object):
    def __init__(self):
        # model parameters
        self.n_release_sites = 20
        self.release_probability = 0.1
        self.mini_amplitude = 50e-6
        self.mini_amplitude_stdev = 20e-6
        self.recovery_tau = 100e-3
        self.measurement_stdev = 100e-6
    
    def measure_likelihood(self, times, amplitudes):
        """Return a measure of the likelihood that *times* and *amplitudes* could be generated
        by a synapse with the current dynamic parameters.
        """
        # state parameters:
        # available_vesicles is a float as a means of avoiding the need to model stochastic vesicle docking;
        # we just assume that recovery is a continuous process. 
        available_vesicles = 10
        
        sample_likelihood = []
        sample_available_vesicles = []
        last_t = times[0]
        for i in range(len(times)):
            t = times[i]
            amplitude = amplitudes[i]
            
            # measure likelihood of seeing this amplitude
            likelihood = self._likelihood([amplitude], available_vesicles)
            sample_likelihood.append(likelihood[0])
            
            # update state
            dt = t - last_t
            last_t = t
            
            # recover vesicles
            recovery = np.exp(-dt / self.recovery_tau)
            available_vesicles = available_vesicles * recovery + self.n_release_sites * (1.0 - recovery)
            
            # release vesicles
            available_vesicles -= amplitude / self.mini_amplitude
            
            sample_available_vesicles.append(available_vesicles)
            
        self.sample_available_vesicles = np.array(sample_available_vesicles)
        
        return np.array(sample_likelihood)
            
    def _likelihood(self, amplitudes, available_vesicles):
        """Estimate the probability density of seeing a particular *amplitude*
        given a number of *available_vesicles*.
        """
        available_vesicles = int(max(0, available_vesicles))
        
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


if __name__ == '__main__':
    
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

    model = StochasticReleaseModel()
    model.mini_amplitude = first_pulse_amps.mean() / (model.n_release_sites * model.release_probability)
    model.mini_amplitude_stdev = model.mini_amplitude / 3.
    
    raw_bg_events = get_baseline_amps(session, pair, clamp_mode='ic')
    mask = event_qc(raw_bg_events)
    bg_events = raw_bg_events[mask]

    model.measurement_stdev = bg_events[amplitude_field].std()
    

    # 3. For each point in the parameter space, simulate a synapse and estimate the joint probability of the set of measured amplitudes
    
    likelihood = model.measure_likelihood(spike_times, amplitudes)


    # 4. Visualize / characterize mapped parameter space. Somehow.
    
    # color events by likelihood
    cmap = pg.ColorMap([0, 1.0], [(0, 0, 0), (255, 0, 0)])
    err = 1.0 / likelihood
    err_norm = 0.5 + (err - err.mean()) / err.std()
    err_colors = cmap.map(err_norm)
    brushes = [pg.mkBrush(c) for c in err_colors]

    # log spike intervals to make visualization a little easier
    compressed_spike_times = np.empty(len(spike_times))
    compressed_spike_times[0] = 0.0
    np.cumsum(np.diff(spike_times)**0.25, out=compressed_spike_times[1:])
    
    win = pg.GraphicsLayoutWidget()
    win.show()
    plt1 = win.addPlot(0, 0, title="likelihood vs compressed time")
    plt1.plot(compressed_spike_times, likelihood, pen=None, symbol='o', symbolBrush=brushes)
    plt2 = win.addPlot(1, 0, title="deconvolved amplitude vs compressed time")
    plt2.plot(compressed_spike_times, amplitudes, pen=None, symbol='o', symbolBrush=brushes)
    plt2.setXLink(plt1)
    plt3 = win.addPlot(2, 0, title="available_vesicles vs compressed time")
    plt3.plot(compressed_spike_times, model.sample_available_vesicles, pen=None, symbol='o', symbolBrush=brushes)
    plt3.setXLink(plt1)



