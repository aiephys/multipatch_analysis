"""Analyses for fitting the average response from a pair given an optional response latency.
Spike time is set to 0ms
"""

import sys
from collections import OrderedDict
import numpy as np
import pyqtgraph as pg
from neuroanalysis.data import TSeries, TSeriesList
from neuroanalysis.baseline import float_mode
from neuroanalysis.fitting import Psp, StackedPsp
from aisynphys.database import default_db as db
import aisynphys.data.data_notes_db as notes_db
from aisynphys.qc import spike_qc
from aisynphys.fitting import fit_avg_pulse_response


def get_pair_avg_fits(pair, session, notes_session=None, ui=None):
    """Return PSP fits to averaged responses for this pair.
    
    Operations are:
    - query all pulse responses for this pair
    - sort responses by clamp mode and holding potential
    - generate average response for each mode/holding combination
    - fit averages to PSP curve
    
    Returns
    -------
    results : dict
        {(mode, holding): {
            'traces': , 
            'average', 
            'fit_params',
            'initial_latency',
            'fit_qc_pass',
            'expected_fit_params',
            'avg_baseline_noise',
            }, 
        }
    
    """
    prof = pg.debug.Profiler(disabled=True, delayed=False)
    prof(str(pair))
    results = {}
    
    # query and sort pulse responses
    records = response_query(session=session, pair=pair).all()
    prof('query prs')
    pulse_responses = [rec[0] for rec in records]
    sorted_responses = sort_responses(pulse_responses)
    prof('sort prs')

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
            latency_window = (0.5e-3, 8e-3)
        else:
            notes = notes_rec.notes
            if notes.get('fit_parameters') is None:
                init_latency = None
                latency_window = (0.5e-3, 8e-3)
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

        # load up expected fit results and compare to manually-verified
        # results
        if notes is None:
            qc_pass = False
            reasons = ['no data notes entry']
            expected_fit_params = None
        elif 'fit_pass' not in notes:
            print("=============================")
            print(pair)
            print(pair.has_synapse)
            print(notes)
            raise Exception("WTF")
        elif notes['fit_pass'][clamp_mode][str(holding)] is not True:
            qc_pass = False
            reasons = ['data notes fit failed qc']
            expected_fit_params = None
        else:
            expected_fit_params = notes['fit_parameters']['fit'][clamp_mode][str(holding)]
            qc_pass, reasons = check_fit_qc_pass(fit_result, expected_fit_params, clamp_mode)

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

    print('; '.join(failures))
    return len(failures) == 0, failures
