import sys
import numpy as np
import pyqtgraph as pg
import scipy.stats as stats
from multipatch_analysis.database import database as db
from multipatch_analysis.pulse_response_strength import PulseResponseStrength
from multipatch_analysis.connection_strength import ConnectionStrength, get_amps
from neuroanalysis.synaptic_release import ReleaseModel


class StochasticReleaseModel(object):
    def __init__(self):
        # model parameters
        self.n_release_sites = 20
        self.release_probability = 0.1
        self.mini_amplitude = 50e-6
        self.mini_amplitude_stdev = 20e-6
        # self.tau_recovery = 100e-3
        self.measurement_stdev = 100e-6
    
    def measure_likelihood(self, times, amplitudes):
        """Return a measure of the likelihood that *times* and *amplitudes* could be generated
        by a synapse with the current dynamic parameters.
        """
        # state parameters:
        available_vesicles = 10
        
        sample_likelihood = []
        for i in range(len(times)):
            t = times[i]
            amplitude = amplitudes[i]
            
            likelihood = self._likelihood(amplitude, available_vesicles)
            sample_likelihood.append(likelihood)
            
        return sample_likelihood
            
    def _likelihood(self, amplitude, available_vesicles):
        """Estimate the probability density of seeing a particular *amplitude*
        given a number of *available_vesicles*.
        """
        likelihood = 0
        release_prob = stats.binom(available_vesicles, self.release_probability)
        for n_vesicles in range(available_vesicles):
            # probability of releasing n_vesicles given available_vesicles and release_probability
            p_n = release_prob.pmf(n_vesicles)
            
            # expected amplitude for n_vesicles
            amp_mean = n_vesicles * self.mini_amplitude
            
            # distribution of amplitudes expected for n_vesicles
            amp_stdev = (self.mini_amplitude_stdev**2 + self.measurement_stdev**2) ** 0.5
            amp_prob = stats.norm(amp_mean, amp_stdev).pdf(amplitude)
            
            # increment likelihood of seeing this amplitude
            likelihood += p_n * amp_prob
        
        return likelihood


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

    events = get_amps(session, pair, clamp_mode='ic')
    mask = events['ex_qc_pass'] == True
    # need more stringent qc for dynamics:
    mask &= np.abs(events['baseline_current']) < 500e-12
    events = events[mask]

    rec_times = 1e-9 * (events['rec_start_time'].astype(float) - float(events['rec_start_time'][0]))
    spike_times = events['max_dvdt_time'] + rec_times

    # log spike intervals to make visualization a little easier
    compressed_spike_times = np.empty(len(spike_times))
    compressed_spike_times[0] = 0.0
    np.cumsum(np.diff(spike_times)**0.25, out=compressed_spike_times[1:])
    pg.plot(compressed_spike_times, events['pos_dec_amp'], pen=None, symbol='o', title="deconvolved amplitude vs compressed time")


    # 2. Initialize model parameters:
    #    - release model with depression, facilitation
    #    - number of synapses, distribution of per-vesicle amplitudes estimated from first pulse CV
    #    - measured distribution of background noise
    #    - parameter space to be searched

    first_pulse_amps = events[events['pulse_number'] == 1]
    first_pulse_stdev = first_pulse_amps['pos_dec_amp'].std()


            
    model = StochasticReleaseModel()


    # 3. For each point in the parameter space, simulate a synapse and estimate the joint probability of the set of measured amplitudes
    likelihood = model.measure_likelihood(spike_times, 

    # 4. Visualize / characterize mapped parameter space. Somehow.


