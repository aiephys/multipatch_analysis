import numpy as np
import warnings, sys

from pyqtgraph.debug import Profiler
from neuroanalysis.data import TSeriesList
from neuroanalysis.fitting import StackedPsp, Psp, fit_psp
from aisynphys.data import PulseResponseList


def fit_avg_pulse_response(pulse_response_list, latency_window, sign, init_params=None, ui=None):
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
    prof = Profiler(disabled=True, delayed=False)
    pair = pulse_response_list[0].pair
    clamp_mode = pulse_response_list[0].recording.patch_clamp_recording.clamp_mode

    # make a list of spike-aligned postsynaptic tseries
    tsl = PulseResponseList(pulse_response_list).post_tseries(align='spike', bsub=True)
    prof('make tseries list')
    
    # average all together
    average = tsl.mean()
    prof('average')
        
    # start with even weighting
    weight = np.ones(len(average))
    
    # boost weight around PSP onset
    onset_start_idx = average.index_at(latency_window[0])
    onset_stop_idx = average.index_at(latency_window[1] + 4e-3) 
    weight[onset_start_idx:onset_stop_idx] = 3.0
    
    # decide whether to mask out crosstalk artifact
    pre_id = int(pair.pre_cell.electrode.ext_id)
    post_id = int(pair.post_cell.electrode.ext_id)
    if abs(pre_id - post_id) < 3:
        # nearby electrodes; mask out crosstalk
        pass
    prof('weights')

    fit = fit_psp(average, search_window=latency_window, clamp_mode=clamp_mode, sign=sign, baseline_like_psp=True, init_params=init_params, fit_kws={'weights': weight})
    prof('fit')
    
    return fit, average
