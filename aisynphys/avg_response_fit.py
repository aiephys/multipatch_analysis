# coding: utf8
import numpy as np
from pyqtgraph.debug import Profiler
import pyqtgraph as pg
from neuroanalysis.fitting import fit_psp
import aisynphys.data.data_notes_db as notes_db
from aisynphys.database import default_db as db
from aisynphys.data import PulseResponseList


def get_pair_avg_fits(pair, session, notes_session=None, ui=None, max_ind_freq=50):
    """Return PSP fits to averaged responses for this pair.

    Fits are performed against average PSPs in 4 different categories: 
    IC -70mV, IC -55mV, VC -70mV, and VC -55mV. All PSPs in these categories are averaged together
    regardless of their position in a pulse train, so we expect the amplitudes of these averages to
    be affected by any short-term depression/facilitation present at the synapse. As such, these fits
    are not ideal for measuring the amplitude of the synapse; however, they do provide good estimates
    of rise time and decay tau.
    
    Operations are:
    
    - Query all pulse responses for this pair, where the pulse train frequency was 
      no faster than max_ind_freq
    - Sort responses by clamp mode and holding potential, with the latter in two bins: -80 to -60 mV and -60 to -45 mV.
      Responses are further separated into qc pass/fail for each bin. QC pass criteria:
        - PR must have exactly one presynaptic spike with detectable latency
        - Either PR.ex_qc_pass or .in_qc_pass must be True, depending on clamp mode / holding
    - Generate average response for qc-passed pulses responses in each mode/holding combination
    - Fit averages to PSP curve. If the latency was manually annotated for this synapse, then the curve
      fit will have its latency constrained within ±100 μs.
    - Compare to manually verified fit parameters; if these are not a close match OR if the 
      manual fits were already failed, then *fit_qc_pass* will be False.
      
    
    Returns
    -------
    results : dict
        {(mode, holding): {
            'responses': ..,
            'average': ..,
            'initial_latency': ..,
            'fit_result': ..,
            'fit_qc_pass': ..,
            'fit_qc_pass_reasons': ..,
            'expected_fit_params': ..,
            'expected_fit_pass': ..,
            'avg_baseline_noise': ..,
        }, ...}
    
    """
    prof = pg.debug.Profiler(disabled=True, delayed=False)
    prof(str(pair))
    results = {}
    
    # query and sort pulse responses with induction frequency 50Hz or slower
    records = response_query(session=session, pair=pair, max_ind_freq=max_ind_freq).all()
    prof('query prs')
    pulse_responses = [rec[0] for rec in records]

    # sort into clamp mode / holding bins
    sorted_responses = sort_responses(pulse_responses)
    prof('sort prs')

    # load expected PSP curve fit parameters from notes DB
    notes_rec = notes_db.get_pair_notes_record(pair.experiment.ext_id, pair.pre_cell.ext_id, pair.post_cell.ext_id, session=notes_session)
    prof('get pair notes')

    if ui is not None:
        ui.show_pulse_responses(sorted_responses)
        ui.show_data_notes(notes_rec)
        prof('update ui')

    for (clamp_mode, holding), responses in sorted_responses.items():
        if len(responses['qc_pass']) == 0:
            results[clamp_mode, holding] = None
            continue
            
        if notes_rec is None:
            notes = None
            sign = 0
            init_latency = None
            latency_window = (0.5e-3, 10e-3)
        else:
            notes = notes_rec.notes
            if notes.get('fit_parameters') is None:
                init_latency = None
                latency_window = (0.5e-3, 10e-3)
            else:
                init_latency = notes['fit_parameters']['initial'][clamp_mode][str(holding)]['xoffset']
                latency_window = (init_latency - 100e-6, init_latency + 100e-6)
            
            # Expected response sign depends on synapse type, clamp mode, and holding:
            sign = 0
            if notes['synapse_type'] == 'ex':
                sign = -1 if clamp_mode == 'vc' else 1
            elif notes['synapse_type'] == 'in' and holding == -55:
                sign = 1 if clamp_mode == 'vc' else -1

        prof('prepare %s %s' % (clamp_mode, holding))
        fit_result, avg_response = fit_avg_pulse_response(responses['qc_pass'], latency_window, sign)
        prof('fit avg')

        # measure baseline noise
        avg_baseline_noise = avg_response.time_slice(avg_response.t0, avg_response.t0+7e-3).data.std()

        # compare to manually-verified results
        if notes is None:
            qc_pass = False
            reasons = ['no data notes entry']
            expected_fit_params = None
            expected_fit_pass = None
        elif notes['fit_pass'][clamp_mode][str(holding)] is not True:
            qc_pass = False
            reasons = ['data notes fit failed qc']
            expected_fit_params = None
            expected_fit_pass = False
        else:
            expected_fit_params = notes['fit_parameters']['fit'][clamp_mode][str(holding)]
            expected_fit_pass = True
            qc_pass, reasons = check_fit_qc_pass(fit_result, expected_fit_params, clamp_mode)
            if not qc_pass:
                print("%s %s %s: %s" % (str(pair), clamp_mode, holding,  '; '.join(reasons)))

        if ui is not None:
            ui.show_fit_results(clamp_mode, holding, fit_result, avg_response, qc_pass)

        results[clamp_mode, holding] = {
            'responses': responses,
            'average': avg_response,
            'initial_latency': init_latency,
            'fit_result': fit_result,
            'fit_qc_pass': qc_pass,
            'fit_qc_pass_reasons': reasons,
            'expected_fit_params': expected_fit_params,
            'expected_fit_pass': expected_fit_pass,
            'avg_baseline_noise': avg_baseline_noise,
        }

    return results


def response_query(session, pair, max_ind_freq=50):
    """Query pulse responses appropriate for generating nice average PSP/PSC shapes.
    
    - Only select from multipatch probes with induction frequencies <= 50Hz
    """
    q = session.query(db.PulseResponse, db.PatchClampRecording, db.StimPulse)
    
    q = q.join(db.StimPulse, db.PulseResponse.stim_pulse)
    q = q.join(db.StimSpike, db.StimSpike.stim_pulse_id==db.StimPulse.id)
    q = q.join(db.Recording, db.PulseResponse.recording)
    q = q.join(db.PatchClampRecording)
    q = q.join(db.MultiPatchProbe)
    
    q = q.filter(db.PulseResponse.pair_id == pair.id)
    q = q.filter(db.MultiPatchProbe.induction_frequency <= max_ind_freq)
    
    return q


def sort_responses(pulse_responses):
    """Sort a list of pulse responses by clamp mode and holding potential into 4 categories: 
    (ic -70), (ic -55), (vc -70), (vc -55). Each category contains pulse responses split into
    lists for qc-pass and qc-fail.

    QC pass for this function requires that the pulse response has True for either pr.in_qc_pass
    or pr.ex_qc_pass, depending on which category the PR has been sorted into. We _also_ require
    at this stage that the PR has exactly one presynaptic spike with a detectable onset time.
    """
    ex_limits = [-80e-3, -60e-3]
    in_limits = [-60e-3, -45e-3]
    
    sorted_responses = {
        ('ic', -70): {'qc_pass': [], 'qc_fail': []},
        ('ic', -55): {'qc_pass': [], 'qc_fail': []},
        ('vc', -70): {'qc_pass': [], 'qc_fail': []},
        ('vc', -55): {'qc_pass': [], 'qc_fail': []},
    }
    qc = {False: 'qc_fail', True: 'qc_pass'}
    
    for pr in pulse_responses: 
        post_rec = pr.recording
        clamp_mode = post_rec.patch_clamp_recording.clamp_mode
        holding = post_rec.patch_clamp_recording.baseline_potential

        if in_limits[0] <= holding < in_limits[1]:
            qc_pass = qc[pr.in_qc_pass and pr.stim_pulse.n_spikes == 1 and pr.stim_pulse.first_spike_time is not None]
            sorted_responses[clamp_mode, -55][qc_pass].append(pr)
        elif ex_limits[0] <= holding < ex_limits[1]:
            qc_pass = qc[pr.ex_qc_pass and pr.stim_pulse.n_spikes == 1 and pr.stim_pulse.first_spike_time is not None]
            sorted_responses[clamp_mode, -70][qc_pass].append(pr)
    
    return sorted_responses


def check_fit_qc_pass(fit_result, expected_params, clamp_mode):
    """Return bool indicating whether a PSP fit result matches a previously accepted result, as well as
    a list of strings describing the reasons, if any.
    """
    failures = []
    fit_params = fit_result.best_values

    # decide on relative and absolute thresholds to use for comparing each parameter
    abs_amp_threshold = 40e-6 if clamp_mode == 'ic' else 1e-12
    if fit_result.nrmse() < expected_params['nrmse']:
        # if the new fit improves NRMSE, then we can be a little more forgiving about not matching the previous fit.
        thresholds = {'amp': (0.3, abs_amp_threshold*1.5), 'rise_time': (1.0, 1e-3), 'decay_tau': (3.0, 5e-3), 'xoffset': (1.0, 150e-6)}
    else:
        thresholds = {'amp': (0.2, abs_amp_threshold), 'rise_time': (0.5, 1e-3), 'decay_tau': (2.0, 5e-3), 'xoffset': (1.0, 150e-6)}

    # compare parameters
    for k, (error_threshold, abs_threshold) in thresholds.items():
        v1 = fit_params[k]
        v2 = expected_params[k]
        if v2 == 0:
            continue
        error = abs(v1-v2) / v2
        # We expect large relative errors when the values are small relative to noise,
        # and large absolute errors otherwise.
        if (error > error_threshold) and (abs(v1 - v2) > abs_threshold):
            failures.append('%s error too large (%s != %s)' % (k, v1, v2))

    return len(failures) == 0, failures


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
    
    if len(tsl) == 0:
        return None, None
    
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
