import numpy as np
import warnings, sys

from pyqtgraph.debug import Profiler
from neuroanalysis.data import TSeriesList
from neuroanalysis.fitting import StackedPsp, Psp, fit_psp



def fit_avg_pulse_response(pulse_response_list, latency_window, sign, ui=None):
    """Generate PSP fit parameters for a list of pulse responses, possibly correcting
    for crosstalk artifacts and gap junctional current during the presynaptic stimulus.
    
    Parameters
    ----------
    pulse_response_list : list
        A list of PulseResponse instances to be time-aligned, averaged, and fit.
    latency_window : (float, float)
        Beginning and end times of a window over which to search for a synaptic response,
        relative to the spike time.
    sign : int
        +1, -1, or 0 indicating the expected sign of the response (see neuroanalysis.fitting.fit_psp)
    
    Returns
    -------
    fit : lmfit ModelResult
        The resulting PSP fit
    average : TSeries
        The averaged pulse response data
    
    """
    prof = Profiler(disabled=False, delayed=False)
    pair = pulse_response_list[0].pair
    clamp_mode = pulse_response_list[0].recording.patch_clamp_recording.clamp_mode

    # make a list of spike-aligned postsynaptic tseries    
    tsl = []
    for pr in pulse_response_list:
        spike_t = pr.stim_pulse.first_spike_time
        if spike_t is None:
            continue
        post_ts = pr.post_tseries
        ts = post_ts.copy(t0=post_ts.t0-spike_t)
        tsl.append(ts)
    tsl = TSeriesList(tsl)
    prof('make tseries list')
    
    # average all together
    average = tsl.mean()
    prof('average')
        
    # start with even weighting
    weight = np.ones(len(average))
    
    # boost weight around PSP onset
    onset_start_idx = average.index_at(latency_window[0])
    onset_stop_idx = average.index_at(latency_window[1] + 2e-3) 
    weight[onset_start_idx:onset_stop_idx] = 3.0
    
    # decide whether to mask out crosstalk artifact
    pre_id = int(pr.pair.pre_cell.electrode.ext_id)
    post_id = int(pr.pair.post_cell.electrode.ext_id)
    if abs(pre_id - post_id) < 3:
        # nearby electrodes; mask out crosstalk
        pass
    prof('weights')

    fit = fit_psp(average, search_window=latency_window, clamp_mode=clamp_mode, sign=sign, fit_kws={'weights': weight})
    prof('fit')
    
    return fit, average
    

def fit_avg_response(traces, mode, latency, sign):
        output_fit_parameters = {}

        if latency is None:
            x_offset = 1e-3
            x_offset_win = [-1e-3, 6e-3]
        else:
            x_offset = latency
            x_offset_win = [-0.1e-3, 0.1e-3]
        
        if len(traces) == 0:
            return output_fit_parameters, x_offset, None
        
        grand_trace = TSeriesList(traces).mean()

        weight = np.ones(len(grand_trace.data))*10.  #set everything to ten initially
        weight[int(1e-3/db.default_sample_rate):int(3e-3/db.default_sample_rate)] = 30.  #area around steep PSP rise 
        ic_weight = weight
        ic_weight[0:int(1e-3/db.default_sample_rate)] = 0.   #area around stim artifact

        mode_params = {
            'vc': {
                'stacked': False,
                'initial_rise': 1e-3,
                'rise_bounds': [0.1e-3, 6e-3],
                'weight': weight
            },
            'ic': {
                'stacked': True,
                'initial_rise': 5e-3,
                'rise_bounds': [1e-3, 30e-3],
                'weight': ic_weight
            }
        }
        
        stacked = mode_params[mode]['stacked']
        initial_rise = mode_params[mode]['initial_rise']
        rise_bounds = mode_params[mode]['rise_bounds']
        weight = mode_params[mode]['weight']
        
        rise_times = list(initial_rise*2.**np.arange(-2, 3, 1))
        x_win = [x_offset + x_offset_win[0], x_offset + x_offset_win[1]]
        x_range = list(np.linspace(x_win[0], x_win[1], 4))

        try:
            fit = fit_psp(grand_trace, 
                mode=mode, 
                sign=sign,
                xoffset=(x_range, x_win[0], x_win[1]),
                rise_time=(rise_times, rise_bounds[0], rise_bounds[1]),
                stacked=stacked,
                fit_kws={'tol': 0.01, 'maxiter': 50},
                
            )
            for param, val in fit.best_values.items():
                output_fit_parameters[param] = val
            output_fit_parameters['yoffset'] = fit.best_values['yoffset']
            output_fit_parameters['nrmse'] = fit.nrmse()
        except:
            print("Error in PSP fit:")
            sys.excepthook(*sys.exc_info())
            return output_fit_parameters, x_offset, None

        return output_fit_parameters, x_offset, fit.best_fit, grand_trace
