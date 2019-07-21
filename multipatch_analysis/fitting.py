import numpy as np
import warnings

from pyqtgraph.debug import Profiler
from neuroanalysis.fitting import StackedPsp, Psp
from .data import MultiPatchProbe, Analyzer, PulseStimAnalyzer


def test_create_all_fit_param_combos():
    """Tests for changes in create_all_fit_param_combos module.
    """
    input={'decay_tau': (0.05, 0.005, 0.5), 
           'rise_time': (0.005, 0.0005, 0.05), 
           'yoffset': (0, -float('inf'), float('inf')), 
           'rise_power': ([1, 2], 'fixed'), 
           'amp': (0.0002, 0, 0.1), 
           'xoffset': ([0.011, 0.014, 0.015], -float('inf'), float('inf'))}
    output=[{'decay_tau': (0.05, 0.005, 0.5), 'rise_time': (0.005, 0.0005, 0.05), 'yoffset': (0, -float('inf'), float('inf')), 'rise_power': (1, 'fixed'), 'amp': (0.0002, 0, 0.1), 'xoffset': (0.011, -float('inf'), float('inf'))}, 
            {'decay_tau': (0.05, 0.005, 0.5), 'rise_time': (0.005, 0.0005, 0.05), 'yoffset': (0, -float('inf'), float('inf')), 'rise_power': (2, 'fixed'), 'amp': (0.0002, 0, 0.1), 'xoffset': (0.011, -float('inf'), float('inf'))}, 
            {'decay_tau': (0.05, 0.005, 0.5), 'rise_time': (0.005, 0.0005, 0.05), 'yoffset': (0, -float('inf'), float('inf')), 'rise_power': (1, 'fixed'), 'amp': (0.0002, 0, 0.1), 'xoffset': (0.014, -float('inf'), float('inf'))}, 
            {'decay_tau': (0.05, 0.005, 0.5), 'rise_time': (0.005, 0.0005, 0.05), 'yoffset': (0, -float('inf'), float('inf')), 'rise_power': (2, 'fixed'), 'amp': (0.0002, 0, 0.1), 'xoffset': (0.014, -float('inf'), float('inf'))}, 
            {'decay_tau': (0.05, 0.005, 0.5), 'rise_time': (0.005, 0.0005, 0.05), 'yoffset': (0, -float('inf'), float('inf')), 'rise_power': (1, 'fixed'), 'amp': (0.0002, 0, 0.1), 'xoffset': (0.015, -float('inf'), float('inf'))}, 
            {'decay_tau': (0.05, 0.005, 0.5), 'rise_time': (0.005, 0.0005, 0.05), 'yoffset': (0, -float('inf'), float('inf')), 'rise_power': (2, 'fixed'), 'amp': (0.0002, 0, 0.1), 'xoffset': (0.015, -float('inf'), float('inf'))}]

    test_out=create_all_fit_param_combos(input)
    
    assert test_out == output


def create_all_fit_param_combos(base_params):
    '''Convert the parameters fed to fit_psp into a list of all possible parameter 
    dictionaries to be fed to PSP() or stackedPSP() for fitting. 
    
    Parameters 
    ----------
    base_params: dictionary
        Each value in the key:value dictionary pair must be a tuple.
        In general the structure of the tuple is of the form, 
        (initial conditions, lower boundary, higher boundary).
        The initial conditions can be either a number or a list 
        of numbers specifying several initial conditions.  The 
        initial condition may also be fixed by replacing the lower 
        higher boundary combination with 'fixed'.
        Note that technically the value of a key:value pair could 
        be a single string (this ends up being evaluated as an 
        expression by lmfit later on).
    
    Returns
    -------
    param_dict_list: list of dictionaries
        Each dictionary contains parameter inputs for one fit_psp run.
        
    Examples:    
    base_params[amplitude]=(10, 0, 20)
    base_params[amplitude]=(10, 'fixed')
    base_params[amplitude]=([5,10, 20], 0, 20)
    base_params[amplitude]=([5,10, 20], 'fixed')
    
    '''
    # need to create all combinations of the initial conditions
    param_dict_list = [{}] #initialize list
    for key, value in base_params.items():
        if isinstance(value[0], list):
            temp_list=[]
            for init_cond in value[0]: #for each initial condition
                temp=[pdl.copy() for pdl in param_dict_list] #copy each parameter dictionary in the list so they do not point back to the original dictionary
                for t in temp:  #for each dictionary in the list 
                    t[key]=tuple([init_cond] +list(value[1:])) #add the key and single initial condition pair
                temp_list=temp_list+temp
            param_dict_list=list(temp_list) #Note this works because the dictionaries are already made immutable above
        else: 
            for pdl in param_dict_list:
                pdl[key]=value
    
    return param_dict_list


def fit_psp(response, 
            mode='ic', 
            sign='any', #Note this will not be used if *amp* input is specified
            method='leastsq', 
            fit_kws=None, 
            stacked=True,
            rise_time_mult_factor=10., #Note this will not be used if *rise_time* input is specified 
            weight='default',
            amp_ratio='default', 
            # the following are parameters that can be fit 
                amp='default',
                decay_tau='default',
                rise_power='default',
                rise_time='default',
                exp_amp='default',
                xoffset='default', 
                yoffset='default',
            ):
    """Fit psp. function to the equation 
    
    This function make assumptions about where the cross talk happens as traces 
    have been aligned to the pulses and response.t0 has been set to zero 
    
    Parameters
    ----------
    response : neuroanalysis.data.Trace class
        Contains data on trace waveform.
    mode : string
        either 'ic' for current clamp or 'vc' for voltage clamp
    sign : string
        Specifies the sign of the PSP deflection.  Must be '+', '-', or any. If *amp* 
        is specified, value will be irrelevant.
    method : string 
        Method lmfit uses for optimization
    rise_time_mult_factor: float
        Parameter that goes into the default calculation rise time.  
        Note that if an input for *rise_time* is provided this input
        will be irrelevant.
    stacked : True or False
        If true, use the StackedPsp function which assumes there is still left
        over voltage decay from previous events.  If False, use Psp function
        which assumes the region of the waveform before the event is at baseline.
    fit_kws : dictionary
        Additional key words that are fed to lmfit
    exp_amp : string
        function that is fed to lmfit
    The parameters below are fed to the psp function. Each value in the 
        key:value dictionary pair must be a tuple.
        In general the structure of the tuple is of the form, 
        (initial conditions, lower boundary, higher boundary).
        The initial conditions can be either a number or a list 
        of numbers specifying several initial conditions.  The 
        initial condition may also be fixed by replacing the lower 
        higher boundary combination with 'fixed'.    
        Examples:    
            amplitude=(10, 0, 20)
            amplitude=(10, 'fixed')
            amplitude=([5,10, 20], 0, 20)
            amplitude=([5,10, 20], 'fixed') 
        xoffset : scalar
            Horizontal shift between begin (positive shifts to the right)
        yoffset : scalar
            Vertical offset
        rise_time : scalar
            Time from beginning of psp until peak
        decay_tau : scalar
            Decay time constant
        amp : scalar
            The peak value of the psp
        rise_power : scalar
            Exponent for the rising phase; larger values result in a slower activation 
        amp_ratio : scalar 
            if *stacked* this is used to set up the ratio between the 
            residual decay amplitude and the height of the PSP.
    
    Returns
    -------
    fit: lmfit.model.ModelResult
        Best fit
    """           
    prof = Profiler(disabled=False, delayed=False)    
    # extracting these for ease of use
    t=response.time_values
    y=response.data
    dt = response.dt
    
    # set initial conditions depending on whether in voltage or current clamp
    # note that sign of these will automatically be set later on based on the 
    # the *sign* input
    if mode == 'ic':
        amp_init = .2e-3
        amp_max = 100e-3
        rise_time_init = 5e-3
        decay_tau_init = 50e-3
    elif mode == 'vc':
        amp_init = 20e-12
        amp_max = 500e-12
        rise_time_init = 1e-3
        decay_tau_init = 4e-3
    else:
        raise ValueError('mode must be "ic" or "vc"')

    # Set up amplitude initial values and boundaries depending on whether *sign* are positive or negative
    if sign == '-':
        amps = (-amp_init, -amp_max, 0)
    elif sign == '+':
        amps = (amp_init, 0, amp_max)
    elif sign =='any':
        warnings.warn("You are not specifying the predicted sign of your psp.  This may slow down or mess up fitting")
        amps = (0, -amp_max, amp_max)
    else:
        raise ValueError('sign must be "+", "-", or "any"')
        
    # initial condition, lower boundry, upper boundry    
    base_params = {
        'xoffset': (14e-3, -float('inf'), float('inf')),
        'yoffset': (0, -float('inf'), float('inf')),
        'rise_time': (rise_time_init, rise_time_init/rise_time_mult_factor, rise_time_init*rise_time_mult_factor),
        'decay_tau': (decay_tau_init, decay_tau_init/10., decay_tau_init*10.),
        'rise_power': (2, 'fixed'),
        'amp': amps
    }
    
    # specify fitting function and set up conditions
    if not isinstance(stacked, bool):
        raise Exception("Stacked must be True or False")
    if stacked:
        psp = StackedPsp()
        base_params.update({
            #TODO: figure out the bounds on these
            'exp_amp': 'amp * amp_ratio',
            'amp_ratio': (0, -100, 100),
        })  
    else:
        psp = Psp()
    
    # override defaults with input
    for bp in base_params.keys():
        if eval(bp)!='default':
            base_params[bp]=eval(bp)
    
    # set weighting that 
    if weight == 'default': #use default weighting
        # THIS CODE IS DEPENDENT ON THE DATA BEING INPUT IN A CERTAIN WAY THAT IS NOT TESTED
        weight = np.ones(len(y))*10.  #set everything to ten initially
        weight[int(10e-3/dt):int(12e-3/dt)] = 0.   #area around stim artifact
        weight[int(12e-3/dt):int(19e-3/dt)] = 30.  #area around steep PSP rise 
    elif weight is False: #do not weight any part of the stimulus
        weight = np.ones(len(y))
    elif 'weight' in vars():  #works if there is a value specified in weight
        if len(weight) != len(y):
            raise Exception('the weight and array vectors are not the same length') 
    
    # arguement to be passed through to fitting function
    fit_kws={'weights': weight}   

    # convert initial parameters into a list of dictionaries to be consumed by psp.fit()        
    param_dict_list= create_all_fit_param_combos(base_params)
    prof('init')

    # cycle though different parameters sets and chose best one
    best_fit = None
    best_score = None
    for p in param_dict_list:
        fit = psp.fit(y, x=t, params=p, fit_kws=fit_kws, method=method)
        prof('fit: %s' % p)
        err = np.sum(fit.residual**2)  # note: using this because normalized (nrmse) is not necessary to comparing fits within the same data set
        if best_fit is None or err < best_score:
            best_fit = fit
            best_score = err
    fit = best_fit

    # nrmse = fit.nrmse()
    if 'baseline_std' in response.meta:
        fit.snr = abs(fit.best_values['amp']) / response.meta['baseline_std']
        fit.err = fit.nrmse() / response.meta['baseline_std']

    return fit
