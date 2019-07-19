"""Analyses for fitting the average response from a pair given an optional response latency.
Spike time is set to 0ms
"""

from multipatch_analysis.database import default_db as db
import multipatch_analysis.data_notes_db as notes_db
import pyqtgraph as pg
import sys
import numpy as np
from collections import OrderedDict
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.baseline import float_mode
from neuroanalysis.fitting import Psp, StackedPsp
from multipatch_analysis.qc import spike_qc
from multipatch_analysis.fitting import fit_psp


def response_query(session, pair):
    q = session.query(
        db.PulseResponse.id.label('response_id'),
        db.PulseResponse.data,
        db.PulseResponse.ex_qc_pass,
        db.PulseResponse.data_start_time.label('rec_start'),
        db.StimPulse.data.label('spike'),
        db.StimPulse.data_start_time.label('spike_start'),
        db.StimPulse.n_spikes,
        db.StimPulse.first_spike_time,
        db.StimSpike.max_slope_time.label('spike_time'),
        db.PatchClampRecording.clamp_mode,
        db.PatchClampRecording.baseline_potential,
        db.MultiPatchProbe.induction_frequency.label('ind_freq'),
    )
    q = q.join(db.StimPulse, db.PulseResponse.stim_pulse)
    q = q.join(db.StimSpike, db.StimSpike.stim_pulse_id==db.StimPulse.id)
    q = q.join(db.Recording, db.PulseResponse.recording)
    q = q.join(db.PatchClampRecording)
    q = q.join(db.MultiPatchProbe)
    q = q.filter(db.PulseResponse.pair_id==pair.id)

    return q


def pair_notes_query(session, pair):
    expt_id = '%0.3f' % pair.experiment.acq_timestamp
    pre_cell_id = str(pair.pre_cell.ext_id)
    post_cell_id = str(pair.post_cell.ext_id)
    q = session.query(notes_db.PairNotes)
    q = q.filter(notes_db.PairNotes.expt_id==expt_id)
    q = q.filter(notes_db.PairNotes.pre_cell_id==pre_cell_id)
    q = q.filter(notes_db.PairNotes.post_cell_id==post_cell_id)

    return q


def sort_responses(pulse_responses):
    ex_limits = [-80e-3, -63e-3]
    in_limits = [-62e-3, -45e-3]
    qc = {False: 'qc_fail', True: 'qc_pass'}
    traces = OrderedDict([
        ('vc', {'-55': {'qc_pass': [], 'qc_fail': []}, '-70': {'qc_pass': [], 'qc_fail': []}}), 
        ('ic', {'-55': {'qc_pass': [], 'qc_fail': []}, '-70': {'qc_pass': [], 'qc_fail': []}}),
    ])
    spikes = OrderedDict([
        ('vc', {'-55': {'qc_pass': [], 'qc_fail': []}, '-70': {'qc_pass': [], 'qc_fail': []}}), 
        ('ic', {'-55': {'qc_pass': [], 'qc_fail': []}, '-70': {'qc_pass': [], 'qc_fail': []}}),
    ])

    for rec in pulse_responses: 
        if rec.ind_freq not in [10, 20, 50]:
            continue
        data_trace = Trace(data=rec.data, sample_rate=db.default_sample_rate)
        spike_trace = Trace(data=rec.spike, sample_rate=db.default_sample_rate)
        n_spikes = rec.n_spikes
        spike_time = rec.spike_time
        if spike_time is None:
            continue 
        clamp = rec.clamp_mode
        holding = rec.baseline_potential
        data_trace.t0 = rec.rec_start - spike_time
        spike_trace.t0 = rec.spike_start - spike_time
        baseline_trace = data_trace.time_slice(-10e-3, -4e-3)
        baseline = float_mode(baseline_trace.data)
        bsub_data_trace = data_trace - baseline
        spike_qc_pass, trace_qc_pass = spike_qc(n_spikes, rec.ex_qc_pass)

        if in_limits[0] < holding < in_limits[1]:
            traces[clamp]['-55'][qc[trace_qc_pass]].append(bsub_data_trace)
            spikes[clamp]['-55'][qc[spike_qc_pass]].append(spike_trace)
        elif ex_limits[0] < holding < ex_limits[1]:
            traces[clamp]['-70'][qc[trace_qc_pass]].append(bsub_data_trace)
            spikes[clamp]['-70'][qc[spike_qc_pass]].append(spike_trace)
    return traces, spikes


def fit_avg_response(traces, mode, holding, latency, sign):
        output_fit_parameters = {}

        if latency is None:
            x_offset = 1e-3
            x_offset_win = [-1e-3, 6e-3]
        else:
            x_offset = latency
            x_offset_win = [-0.1e-3, 0.1e-3]
        
        if len(traces[mode][holding]['qc_pass']) == 0:
            return output_fit_parameters, x_offset, None
        
        grand_trace = TraceList(traces[mode][holding]['qc_pass']).mean()

        weight = np.ones(len(grand_trace.data))*10.  #set everything to ten initially
        weight[int(1e-3/db.default_sample_rate):int(3e-3/db.default_sample_rate)] = 30.  #area around steep PSP rise 
        ic_weight = weight
        ic_weight[0:int(1e-3/db.default_sample_rate)] = 0.   #area around stim artifact

        mode_params = {'vc': {
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
        
        rise_times = list(initial_rise*2.**np.arange(-2, 3, 0.5))
        x_win = [x_offset + x_offset_win[0], x_offset + x_offset_win[1]]
        x_range = list(np.linspace(x_win[0], x_win[1], 7))

        try:
            fit = fit_psp(grand_trace, 
                mode=mode, 
                sign=sign,
                xoffset=(x_range, x_win[0], x_win[1]),
                rise_time=(rise_times, rise_bounds[0], rise_bounds[1]),
                stacked=stacked,
                )
            for param, val in fit.best_values.items():
                output_fit_parameters[param] = val
            output_fit_parameters['yoffset'] = fit.best_values['yoffset']
            output_fit_parameters['nrmse'] = fit.nrmse()
        except:
            print("Error in PSP fit:")
            sys.excepthook(*sys.exc_info())
            return output_fit_parameters, x_offset, None

        return output_fit_parameters, x_offset, fit.best_fit
