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
from multipatch_analysis.database import default_db as db
import multipatch_analysis.data.data_notes_db as notes_db
from multipatch_analysis.qc import spike_qc
from multipatch_analysis.fitting import fit_avg_pulse_response


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
            }, 
        }
    
    """
    prof = pg.debug.Profiler(disabled=False, delayed=False)
    prof(str(pair))
    results = {}
    
    # query and sort pulse responses
    records = response_query(session=session, pair=pair).all()
    prof('query prs')
    pulse_responses = [rec[0] for rec in records]
    sorted_responses = sort_responses(pulse_responses)
    prof('sort prs')

    notes = notes_db.get_pair_notes(pair.experiment.ext_id, pair.pre_cell.ext_id, pair.post_cell.ext_id, session=notes_session)
    prof('get pair notes')

    if ui is not None:
        ui.show_pulse_responses(sorted_responses)
        ui.show_data_notes(notes)
        prof('update ui')

    for (clamp_mode, holding), responses in sorted_responses.items():
        if len(responses['qc_pass']) == 0:
            results[clamp_mode, holding] = None
            continue
            
        if notes is None:
            init_latency = None
            latency_window = (0.5e-3, 8e-3)
            sign = 0
        else:
            init_latency = notes['fit_parameters']['initial'][clamp_mode][str(holding)]['xoffset']
            latency_window = (init_latency - 100e-6, init_latency + 100e-6)
            
            # Expected response sign depends on synapse type, clamp mode, and holding:
            if notes['synapse_type'] == 'ex':
                sign = -1 if clamp_mode == 'vc' else 1
            elif notes['synapse_type'] == 'in':
                if holding == -70:
                    sign = 0
                else:
                    sign = 1 if clamp_mode == 'vc' else -1

        prof('prepare %s %s' % (clamp_mode, holding))
        fit_result, avg_response = fit_avg_pulse_response(responses['qc_pass'], latency_window, sign)
        prof('fit avg')

        if ui is not None:
            ui.show_fit_results(clamp_mode, holding, fit_result, avg_response)

        results[clamp_mode, holding] = {
            'responses': responses,
            'average': avg_response,
            'initial_latency': init_latency,
            'fit_result': fit_result,
            'fit_qc_pass': check_fit_qc_pass(fit_result, notes),
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
    ex_limits = [-80e-3, -63e-3]
    in_limits = [-63e-3, -45e-3]
    
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


def check_fit_qc_pass(fit_result, notes_record):
    if notes_record is None:
        return False
    raise Exception("check notes structure")
    if notes_record['qc_pass'] is not True:
        return False
    expected_params = notes_record['fit_result']
    fit_params = fit_result.best_values
    for k,v1 in fit_params.items():
        v2 = expected_params[k]
        if abs(v1-v2) / v2 > 0.1:
            return False
    if abs(fit_params['xoffset'] - expected_params['xoffset']) > 100e-6:
        return False

    return True
