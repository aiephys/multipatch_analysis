# coding: utf8
from __future__ import print_function, division
import functools
import numpy as np
import numba
import scipy.stats as stats
import scipy.optimize


class StochasticReleaseModel(object):
    """A model of stochastic synaptic release used for determining optimal model parameters that describe
    empirically measured synaptic responses.
    
    Synaptic strength changes moment to moment based on the prior history of action potentials at 
    the presynaptic terminal. Many models have been published previously that attempt to capture this relationship.
    However, synaptic strength also depends on the synapse's prior history of vesicle release. This model uses
    both spike timing and response amplitude to compare a series of evoked synaptic events against the distribution
    of likely amplitudes predicted by the model.
    
    Usually, synaptic response data alone is not sufficient to fully constrain all parameters in a release model.
    This model is intended to be used to search large paremeter spaces to determine the subspace of parameters
    that are consistent with the measured data.
    
    Parameters
    ----------
    params : dict
        A dictionary of parameters specifying the behavior of the model:

        - n_release_sites (int) : Number of synaptic release zones
        - base_release_probability (float) : Resting-state synaptic release probability (0.0-1.0)
        - mini_amplitude (float) : Mean PSP amplitude evoked by a single vesicle release
        - mini_amplitude_cv (float) : Coefficient of variation of PSP amplitude evoked from single vesicle releases
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
        'mini_amplitude_cv',
        'vesicle_recovery_tau',
        'facilitation_amount',
        'facilitation_recovery_tau',
        'measurement_stdev',
    ]
        
    _optimization_plot = None
    
    def __init__(self, params):
        for k in params:
            if k not in self.param_names:
                raise ValueError("Unknown parameter name %r" % k)
        self.params = params

        # How long to wait after a NaN event before the model begins accumulating likelihood values again
        self.missing_event_penalty = 0.0
        

    def optimize_mini_amplitude(self, spike_times, amplitudes, show=False):
        """Given a set of spike times and amplitudes, optimize the mini_amplitude parameter
        to produce the highest likelihood model.
        
        Returns the output of measure_likelihood for the best model.
        """
        params = self.params.copy()

        
        init_amp = estimate_mini_amplitude(amplitudes, params)
        params['mini_amplitude'] = init_amp
        init_result = self.measure_likelihood(spike_times, amplitudes, params)
        mean_amp = np.nanmean(amplitudes)
        if show:
            print("========== Optimize mini amplitude ==============")
            for k,v in params.items():
                print("   %s: %s" % (k, v))
            print("   initial guess:", init_amp)
            print("   mean amp:", mean_amp)
        ratio = mean_amp / np.nanmean(init_result['result']['expected_amplitude'])
        init_amp *= ratio
        if show:
            print("   corrected amp 1:", init_amp)
        
        # in some cases, init_amp is way too large; force these back down:
        init_amp = min(init_amp, mean_amp) if mean_amp > 0 else max(init_amp, mean_amp)
        if show:
            print("   corrected amp 2:", init_amp)
             
        # self.params['mini_amplitude'] = params['mini_amplitude']
        result = self.optimize(spike_times, amplitudes, optimize={'mini_amplitude': (init_amp, init_amp*0.01, init_amp*100)}, show=show)
        if show:
            print("   optimized amp:", result['params']['mini_amplitude'])
        result['optimization_info'] = {'init_amp': init_amp / ratio, 'ratio': ratio, 'corrected_amp': init_amp, 'init_likelihood': init_result['likelihood']}
        return result
    
    def optimize(self, spike_times, amplitudes, optimize, show=False):
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
        best_result['optimization_path'] = {
            'mini_amplitude': [k[0] for k in results.keys()], 
            'likelihood': [v['likelihood'] for v in results.values()]
        }
        
        # update attributes with best result
        for i,k in enumerate(optimize.keys()):
            self.params[k] = best.x[i]
                    
        return best_result
    
    def measure_likelihood(self, spike_times, amplitudes, params=None, show=False):
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
        likelihood = np.exp(np.nanmean(np.log(result['likelihood'] + 0.1)))
        
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
        return release_likelihood(amplitudes, available_vesicles, state['release_probability'], params['mini_amplitude'], params['mini_amplitude_cv'], params['measurement_stdev'])


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
                    mini_amplitude_cv,
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
                likelihood = release_likelihood_scalar(amplitude, av, release_probability, mini_amplitude, mini_amplitude_cv, measurement_stdev)
                # prof('likelihood')
                
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
def release_likelihood(amplitudes, available_vesicles, release_probability, mini_amplitude, mini_amplitude_cv, measurement_stdev):
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
    mini_amplitude_cv : float
        Coefficient of variation of response amplitudes evoked by a single vesicle release
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
       function where µ = nR * mini_amplitude and σ = sqrt((mini_amplitude * mini_amplitude_cv)^2 * nR + measurement_stdev)
    3. The total likelihood is the sum of likelihoods for all possible values of nR.
    """
    return np.array([
        release_likelihood_scalar(amplitude, available_vesicles, release_probability, mini_amplitude, mini_amplitude_cv, measurement_stdev) 
        for amplitude in amplitudes])


@numba.jit(nopython=True)
def release_likelihood_scalar(amplitude, available_vesicles, release_probability, mini_amplitude, mini_amplitude_cv, measurement_stdev):
    """Same as release_likelihood, but optimized for a scalar amplitude argument"""
    n_vesicles = np.arange(available_vesicles + 1)
    
    # probability of releasing n_vesicles given available_vesicles and release_probability
    p_n = binom_pmf_range(available_vesicles, release_probability, available_vesicles + 1)
    
    # expected amplitude for n_vesicles
    amp_mean = n_vesicles * mini_amplitude
    
    # amplitude stdev increases by sqrt(n) with number of released vesicles
    amp_stdev = ((mini_amplitude * mini_amplitude_cv)**2 * n_vesicles + measurement_stdev**2) ** 0.5
    
    # distributions of amplitudes expected for n_vesicles
    amp_prob = p_n * normal_pdf(amp_mean, amp_stdev, amplitude)
    
    # sum all distributions across n_vesicles
    likelihood = amp_prob.sum()
    
    assert likelihood >= 0
    return likelihood


# def release_distribution(available_vesicles, release_probability, mini_amplitude, mini_amplitude_cv, measurement_stdev):
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
#     return amplitudes * sign, release_likelihood(amplitudes, available_vesicles, release_probability, mini_amplitude, mini_amplitude_cv, measurement_stdev) * da

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


def estimate_mini_amplitude(amplitudes, params):
    avg_amplitude = np.nanmean(amplitudes)
    expected = release_expectation_value(params['n_release_sites'], params['base_release_probability'], 1)
    init_amp = avg_amplitude / expected
    
    # takes care of cases where the first peak in the distribution is larger than any measured events;
    # otherwise this case is very difficult to optimize
    while abs(init_amp) > abs(avg_amplitude):
        init_amp /= 2

    return init_amp

