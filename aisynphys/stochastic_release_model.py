# coding: utf8
from __future__ import print_function, division
import functools, pickle, time, os, multiprocessing, traceback
from collections import OrderedDict
import numpy as np
import numba
import scipy.stats as stats
import scipy.optimize


# lets us quickly disable jit for debugging:
def _fake_jit(**kwds):
    return lambda fn: fn
#jit = _fake_jit    
jit = numba.jit


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
        - desensitization_amount (float) : Amount of postsynaptic receptor desensitization to apply (per release site) following each spike
        - desensitization_recovery_tau (float) : Time constant for recovery from desensitization
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
        ('sensitization', float),
    ]

    param_names = [
        'n_release_sites',
        'base_release_probability',
        'mini_amplitude',
        'mini_amplitude_cv',
        'vesicle_recovery_tau',
        'facilitation_amount',
        'facilitation_recovery_tau',
        'desensitization_amount',
        'desensitization_recovery_tau',
        'measurement_stdev',
    ]
        
    def __init__(self, params):
        for k in params:
            if k not in self.param_names:
                raise ValueError("Unknown parameter name %r" % k)
        self.params = params

        # How long to wait after a NaN event before the model begins accumulating likelihood values again
        self.missing_event_penalty = 0.0

    def optimize_mini_amplitude(self, spike_times, amplitudes, event_meta=None, show=False):
        """Given a set of spike times and amplitudes, optimize the mini_amplitude parameter
        to produce the highest likelihood model.
        
        Returns the output of run_model for the best model.
        """
        params = self.params.copy()
        
        init_amp = estimate_mini_amplitude(amplitudes, params)
        params['mini_amplitude'] = init_amp
        init_result = self.run_model(spike_times, amplitudes, params, event_meta=event_meta)
        mean_amp = np.nanmean(amplitudes)
        if show:
            print("========== Optimize mini amplitude ==============")
            for k,v in params.items():
                print("   %s: %s" % (k, v))
            print("   initial guess:", init_amp)
            print("   mean amp:", mean_amp)
        ratio = mean_amp / np.nanmean(init_result.result['expected_amplitude'])
        init_amp *= ratio
        if show:
            print("   corrected amp 1:", init_amp)
        
        # in some cases, init_amp is way too large; force these back down:
        init_amp = min(init_amp, mean_amp) if mean_amp > 0 else max(init_amp, mean_amp)
        if show:
            print("   corrected amp 2:", init_amp)
             
        # self.params['mini_amplitude'] = params['mini_amplitude']
        result = self.optimize(spike_times, amplitudes, optimize={'mini_amplitude': (init_amp, init_amp*0.01, init_amp*100)}, event_meta=event_meta)
        if show:
            print("   optimized amp:", result.optimized_params['mini_amplitude'])
        result.optimization_info = {'init_amp': init_amp / ratio, 'ratio': ratio, 'corrected_amp': init_amp, 'init_likelihood': init_result.likelihood}
        return result
    
    def optimize(self, spike_times, amplitudes, optimize, event_meta=None):
        """Optimize specific parameters to maximize the model likelihood.

        This method updates the attributes for any optimized parameters and returns the
        best result from run_model().
        
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
            res = self.run_model(spike_times, amplitudes, opt_params, event_meta=event_meta)
            results[tuple(x)] = res
            # print(opt_params)
            # print(res['likelihood'])
            return -res.likelihood
        
        best = scipy.optimize.minimize(fn, x0=init, 
            method="Nelder-Mead", options={'fatol': 0.01}  # no bounds, can't set initial step?
            #method='BFGS', bounds=bounds, options={'gtol': 1, 'eps': 10e-6}  # oscillates
            #method='Powell', options={'ftol': 0.01}  # not efficient; can't set initial step
            #method='CG', options={'eps': 10e-6}  # not efficient
            #method='Newton-CG',  # requires jacobian
            #method='L-BFGS-B'  # not efficient
            #method='COBYLA', options={'rhobeg': -100e-6}  # fast but oscillates, misses peak
        )
        
        best_result = results[tuple(best.x.flatten())]

        # take optimized params out of result.params and put them in result.optimized_params instead
        # (to ensure we can re-run in the same way)
        best_result.params = self.params.copy()
        best_result.optimized_params = {k:best.x[i] for i,k in enumerate(optimize.keys())}
        best_result.optimization_init = optimize
        best_result.optimization_result = best
        best_result.optimization_path = {
            'mini_amplitude': [k[0] for k in results.keys()], 
            'likelihood': [v.likelihood for v in results.values()]
        }
        
        # update attributes with best result
        self.params = self.params.copy()
        for i,k in enumerate(optimize.keys()):
            self.params[k] = best.x[i]
                    
        return best_result
    
    def run_model(self, spike_times, amplitudes, params=None, event_meta=None):
        """Run the stochastic release model with a specific set of spike times.
        
        This can be used two different ways: (1) compute a measure of the likelihood that *times* and *amplitudes*
        could be generated by a synapse with the current dynamic parameters, or (2) run the model to predict
        the behavior of a synapse given a set of presynaptic spike times.
        
        Parameters
        ----------
        spike_times : array
            Times (in seconds) of presynaptic spikes in ascending order
        amplitudes : array | str
            If this argument is an array, then it specifies the evoked PSP/PSC amplitudes for each spike 
            listed in *spike_times*. Amplitudes may be NaN to indicate that the event should be ignored 
            (usually because the spike could not be detected or the amplitude could not be measured). 
            Any events within 10 seconds following a NaN will update the model as usual, but will not 
            be included in the likelihood estimate.
            If this argument is a string, then the model will generate simulated amplitudes as it runs.
            The string may be "random" to randomly select amplitudes from the model distribution, or
            "expected" to use the expectation value of the predicted distribution at each spike time.
        params :  dict
            Dictionary of model parameter values to use during this run. By default, parameters
            are taken from self.params.
        event_meta : array
            Extra per-event metadata to be included in the model results
        
        Returns
        -------
        result : StochasticReleaseModelResult
        """
        if params is None:
            params = self.params
        use_expectation = False
        if isinstance(amplitudes, str):
            assert amplitudes in ('expected', 'random'), "amplitudes argument must be ndarray, 'expected', or 'random'"
            use_expectation = amplitudes == 'expected'
            amplitudes = np.array([])
        
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
            use_expectation=use_expectation,
            **params,
        )
        
        # scalar representation of overall likelihood
        likelihood = np.exp(np.nanmean(np.log(result['likelihood'] + 0.1)))
        
        return StochasticReleaseModelResult(
            result=result, 
            pre_spike_state=pre_spike_state,
            post_spike_state=post_spike_state,
            likelihood=likelihood,
            params=params,
            model=self,
            event_meta=event_meta,
        )
            
    def likelihood(self, amplitudes, state, params=None):
        """Estimate the probability density of seeing a particular *amplitude*
        given a number of *available_vesicles*.
        """
        if params is None:
            params = self.params.copy()
        
        available_vesicles = int(np.clip(np.round(state['available_vesicle']), 0, params['n_release_sites']))
        return release_likelihood(amplitudes, available_vesicles, state['release_probability'], params['mini_amplitude'], params['mini_amplitude_cv'], params['measurement_stdev'])

    @staticmethod
    @jit(nopython=True)
    def _run_model( spike_times, 
                    amplitudes,
                    result,
                    pre_spike_state,
                    post_spike_state, 
                    missing_event_penalty,
                    use_expectation,
                    n_release_sites,
                    base_release_probability,
                    mini_amplitude,
                    mini_amplitude_cv,
                    vesicle_recovery_tau,
                    facilitation_amount,
                    facilitation_recovery_tau,
                    desensitization_amount,
                    desensitization_recovery_tau,
                    measurement_stdev,
                    ):

        # initialize state parameters:
        # available_vesicles is a float as a means of avoiding the need to model stochastic vesicle docking;
        # we just assume that recovery is a continuous process. (maybe we should just call this "neurotransmitter"
        # instead)
        
        have_amps = len(amplitudes) > 0
        
        available_vesicle = n_release_sites
        release_probability = base_release_probability
        sensitization = 1.0
        
        previous_t = spike_times[0]
        last_nan_time = -np.inf

        for i,t in enumerate(spike_times):
            if have_amps:
                amplitude = amplitudes[i]
            
            if have_amps and np.isnan(amplitude):
                expected_amplitude = np.nan
                likelihood = np.nan
                last_nan_time = t

            else:
                dt = t - previous_t
                previous_t = t

                # recover vesicles up to the current timepoint
                v_recovery = np.exp(-dt / vesicle_recovery_tau)
                available_vesicle += (n_release_sites - available_vesicle) * (1.0 - v_recovery)

                # apply recovery from facilitation toward baseline release probability
                f_recovery = np.exp(-dt / facilitation_recovery_tau)
                release_probability += (base_release_probability - release_probability) * (1.0 - f_recovery)
                # prof('recover facilitation')

                # recover from desensitization
                s_recovery = np.exp(-dt / desensitization_recovery_tau)
                sensitization += (1.0 - sensitization) * (1.0 - s_recovery)
                
                # predict most likely amplitude for this spike (just for show)
                effective_available_vesicle = max(0, available_vesicle)
                effective_mini_amplitude = mini_amplitude * sensitization
                expected_amplitude = release_expectation_value(
                    effective_available_vesicle,
                    release_probability, 
                    effective_mini_amplitude,
                )
                
                if not have_amps:
                    if use_expectation:
                        # run model forward based on expectation value
                        amplitude = expected_amplitude
                    else:
                        # select a random amplitude from distribution
                        amplitude = release_random_value(
                            effective_available_vesicle,
                            release_probability, 
                            effective_mini_amplitude,
                            mini_amplitude_cv,
                            measurement_stdev
                        )

                pre_available_vesicle = available_vesicle
                pre_release_probability = release_probability
                pre_sensitization = sensitization

                # measure likelihood of seeing this response amplitude
                av = max(0, min(n_release_sites, int(np.round(available_vesicle))))
                likelihood = release_likelihood_scalar(amplitude, av, release_probability, mini_amplitude, mini_amplitude_cv, measurement_stdev)
                # prof('likelihood')
                
                # release vesicles
                # note: we allow available_vesicle to become negative because this helps to ensure
                # that the overall likelihood will be low for such models
                depleted_vesicle = amplitude / mini_amplitude
                available_vesicle -= depleted_vesicle

                assert np.isfinite(available_vesicle)
                
                # apply spike-induced facilitation in release probability
                release_probability += (1.0 - release_probability) * facilitation_amount
                # prof('update state')

                # apply spike-induced desensitization of postsynaptic receptors
                sensitization *= (1.0 - desensitization_amount) ** depleted_vesicle

                # ignore likelihood for this event if it was too close to an unmeasurable response
                if t - last_nan_time < missing_event_penalty:
                    likelihood = np.nan

                # record model state immediately before spike
                pre_spike_state[i]['available_vesicle'] = pre_available_vesicle
                pre_spike_state[i]['release_probability'] = pre_release_probability
                pre_spike_state[i]['sensitization'] = pre_sensitization
                    
                # record model state immediately after spike
                post_spike_state[i]['available_vesicle'] = available_vesicle
                post_spike_state[i]['release_probability'] = release_probability
                post_spike_state[i]['sensitization'] = sensitization
                # prof('record')
            
            if np.isnan(available_vesicle):
                raise Exception("NaNs where they shouldn't be")
            
            # record results
            result[i]['spike_time'] = t
            result[i]['amplitude'] = amplitude
            result[i]['expected_amplitude'] = expected_amplitude
            result[i]['likelihood'] = likelihood


class StochasticReleaseModelResult:
    """Contains the results of StochasticReleaseModel.run_model in several attributes:

    Attributes
    ----------
    result : array
        Array of model results, one record per input spike.
        Contains fields: spike_time, amplitude, expected_amplitude, likelihood
    pre_spike_state : array
        State variable values immediately before each spike
    post_spike_state : array
        State variable values immediately after each spike
    likelihood : float
        The estimated likelihood for this model run
    params : dict
        A dictionary of parameters used to generate the model result
    model : StochasticReleaseModel
        The model instance that generated this result
    event_meta : array
        Per-event metadata, mainly regarding stimulus structure
    optimized_params : dict | None
        Parameters that were optimized by running the model
    """
    pickle_attributes = ['likelihood', 'params', 'optimized_params']

    def __init__(self, result, pre_spike_state, post_spike_state, likelihood, params, model, optimized_params=None, event_meta=None):
        self.result = result
        self.pre_spike_state = pre_spike_state
        self.post_spike_state = post_spike_state
        self.likelihood = likelihood
        self.params = params
        self.model = model
        self.event_meta = event_meta

        # filled in later by optimization routine
        self.optimized_params = optimized_params or {}

    @property
    def all_params(self):
        """A dictionary containing the combination of self.params and self.optimized_params; can be used 
        to re-run the model with the result of the parameter optimization fixed.
        """
        p = self.params.copy()
        p.update(self.optimized_params)
        return p

    def __getstate__(self):
        return {k: getattr(self, k) for k in self.pickle_attributes}

    def __setstate__(self, state):
        for k,v in state.items():
            setattr(self, k, v)


@jit(nopython=True)
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


@jit(nopython=True)
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

@jit(nopython=True)
def release_expectation_value(available_vesicles, release_probability, mini_amplitude):
    """Return the expectation value for the release amplitude distribution defined by the arguments.
    """
    return binom_mean(available_vesicles, release_probability) * mini_amplitude


@jit(nopython=True)
def release_random_value(
    available_vesicle,
    release_probability, 
    mini_amplitude,
    mini_amplitude_cv,
    measurement_stdev):
        n_vesicles = np.random.binomial(n=int(available_vesicle), p=release_probability)
        amp_mean = n_vesicles * mini_amplitude
        amp_stdev = ((mini_amplitude * mini_amplitude_cv)**2 * n_vesicles + measurement_stdev**2) ** 0.5
        return np.random.normal(loc=amp_mean, scale=amp_stdev)


@jit(nopython=True)
def normal_pdf(mu, sigma, x):
    """Probability density function of normal distribution
    
    Same as scipy.stats.norm(mu, sigma).pdf(x)
    """
    return (1.0 / (2 * np.pi * sigma**2))**0.5 * np.exp(- (x-mu)**2 / (2 * sigma**2))

#@functools.lru_cache(maxsize=2**14)
@jit(nopython=True)
def binom_pmf_range(n, p, k):
    """Probability mass function of binomial distribution
    
    Same as scipy.stats.binom(n, p).pmf(arange(k)), but much faster.
    """
    k = np.arange(k)
    bc = np.array([binom_coeff(n,k1) for k1 in k])
    return bc * p**k * (1-p)**(n-k)


_binom_coeff_cache = np.fromfunction(scipy.special.binom, (67, 67)).astype(int)

@jit(nopython=True)
def binom_coeff(n, k):
    """Binomial coefficient: n! / (k! (n-k)!)
    
    Same as scipy.special.binom, but much faster and limited to n < 67.
    """
    # note: one cold imagine writing an optimized binomial coefficient function that
    # is not limited to n < 67, but having to look out for integer overflows slows us down.
    return _binom_coeff_cache[n, k]


@jit(nopython=True)
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


class ParameterSpace(object):
    """Used to generate and store model results over a multidimentional parameter space.

    Parameters
    ----------
    params : dict
        Dictionary of {'parameter_name': array_of_values} describing parameter space to be searched.
        Scalar values are passed to the evaluation function when calling run(), but are not included
        as an axis in the result array.

    """
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
        """Run *func* in parallel over the entire parameter space, storing
        results into self.result.

        If workers==1, then run locally to make debugging easier.
        """
        if workers is None:
            workers = multiprocessing.cpu_count()
        all_inds = list(np.ndindex(self.result.shape))
        all_params = [self[inds] for inds in all_inds] 
        if workers > 1:
            # from pyqtgraph.multiprocess import Parallelize
            # with Parallelize(enumerate(all_inds), results=self.result, progressDialog={'labelText':'synapticulating...', 'nested':True}, workers=workers) as tasker:
            #     for i, inds in tasker:
            #         params = self[inds]
            #         tasker.results[inds] = func(params, **kwds)
            pool = multiprocessing.Pool(workers)
            fn = functools.partial(func, **kwds)
            import pyqtgraph as pg
            with pg.ProgressDialog('synapticulating...', maximum=len(all_inds)) as dlg:
                for i,r in enumerate(pool.imap(fn, all_params, chunksize=100)):
                    dlg += 1
                    if dlg.wasCanceled():
                        pool.terminate()
                        raise Exception("Synapticulation cancelled. No refunds.")
                    self.result[all_inds[i]] = r
        else:
            import pyqtgraph as pg
            with pg.ProgressDialog('synapticulating (serial)...', nested=True, maximum=len(all_inds)) as dlg:
                for inds in all_inds:
                    params = self[inds]
                    self.result[inds] = func(params, **kwds)
                    dlg += 1
                    assert not dlg.wasCanceled()
        
    def __getitem__(self, inds):
        """Return a dict of the parameter values used at a specific index in the parameter space.
        """
        params = self.static_params.copy()
        for i,param in enumerate(self.param_order):
            params[param] = self.params[param][inds[i]]
        return params


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
        db.MultiPatchProbe.induction_frequency,
        db.MultiPatchProbe.recovery_delay,
        db.SyncRec.ext_id.label('sync_rec_ext_id'),
    )
    q = q.join(db.Baseline, db.PulseResponse.baseline)
    q = q.join(db.PulseResponseFit)
    q = q.join(db.StimPulse)
    q = q.join(db.Recording, db.PulseResponse.recording)
    q = q.join(db.SyncRec, db.Recording.sync_rec)
    q = q.join(db.PatchClampRecording)
    q = q.join(db.MultiPatchProbe)

    q = q.filter(db.PulseResponse.pair_id==pair.id)
    q = q.filter(db.PatchClampRecording.clamp_mode=='ic')
    
    q = q.order_by(db.Recording.start_time).order_by(db.StimPulse.onset_time)

    return q


class StochasticModelRunner:
    """Handles loading data for a synapse and executing the model across a parameter space.
    """
    def __init__(self, db, experiment_id, pre_cell_id, post_cell_id, workers=None):
        self.db = db
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
        self.run_model(param_space[(0,) * len(search_params)])

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
        session = self.db.session()

        expt = self.db.experiment_from_ext_id(self.experiment_id, session=session)
        pair = expt.pairs[(self.pre_cell_id, self.post_cell_id)]

        syn_type = pair.synapse.synapse_type
        print("Synapse type:", syn_type)

        # 1. Get a list of all presynaptic spike times and the amplitudes of postsynaptic responses

        events = event_query(pair, self.db, session).dataframe()
        print("loaded %d events" % len(events))

        if len(events) == 0:
            raise Exception("No events found for this synapse.")

        rec_times = (events['rec_start_time'] - events['rec_start_time'].iloc[0]).dt.total_seconds().to_numpy()
        spike_times = events['first_spike_time'].to_numpy() + rec_times
        
        # some metadata to follow the events around--not needed for the model, but useful for 
        # analysis later on.
        event_meta = events[['sync_rec_ext_id', 'pulse_number', 'induction_frequency', 'recovery_delay']]
        
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
        
        return spike_times, amplitudes, bg_amplitudes, event_meta

    @property
    def parameters(self):
        """A structure defining the parameters to search.
        """
        if self._parameters is None:
            self._parameters = self._generate_parameters()
        return self._parameters  

    def _generate_parameters(self):
        spike_times, amplitudes, bg_amplitudes, event_meta = self.synapse_events
        
        search_params = {
            # If mini_amplitude is commented out here, then it will be optimized automatically by the model:
            #'mini_amplitude': np.nanmean(amplitudes) * 1.2**np.arange(-12, 24, 2),

            'n_release_sites': np.array([1, 2, 4, 8, 16, 32, 64]),
            'base_release_probability': np.array([0.00625, 0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]),
            'mini_amplitude_cv': np.array([0.05, 0.1, 0.2, 0.4, 0.8]),
            'measurement_stdev': np.nanstd(bg_amplitudes),

            'vesicle_recovery_tau': np.array([0.0025, 0.01, 0.04, 0.16, 0.64, 2.56]),
            # 'vesicle_recovery_tau': np.array([0.0001, 0.01]),

            'facilitation_amount': np.array([0.0, 0.00625, 0.025, 0.05, 0.1, 0.2, 0.4]),
            'facilitation_recovery_tau': np.array([0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28]),
            # 'facilitation_amount': np.array([0, 0.5]),
            # 'facilitation_recovery_tau': np.array([0.1, 0.5]),

            # 'desensitization_amount': np.array([0.0, 0.00625, 0.025, 0.05, 0.1, 0.2, 0.4]),
            # 'desensitization_recovery_tau': np.array([0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28]),
            # 'desensitization_amount': np.array([0.0, 0.1]),
            # 'desensitization_recovery_tau': np.array([0.01, 0.1]),
            'desensitization_amount': 0,
            'desensitization_recovery_tau': 1,
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

    def run_model(self, params, **kwds):
        model = StochasticReleaseModel(params)
        spike_times, amplitudes, bg, event_meta = self.synapse_events
        if 'mini_amplitude' in params:
            return model.run_model(spike_times, amplitudes, event_meta=event_meta, **kwds)
        else:
            return model.optimize_mini_amplitude(spike_times, amplitudes, event_meta=event_meta, **kwds)


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
        
    def run_model(self, params):
        params = params.copy()
        runner = self.model_runners[params.pop('synapse')]
        return runner.run_model(params)
