# coding: utf8
"""
Analyses that measure the strength of synaptic connections.

"""
from __future__ import print_function, division

import sys, multiprocessing

import numpy as np
import scipy.stats

from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import Trace, TraceList
from neuroanalysis import filter
from neuroanalysis.event_detection import exp_deconvolve
from neuroanalysis.baseline import float_mode
from neuroanalysis.fitting import Psp

from multipatch_analysis.database import database as db
from multipatch_analysis.connection_detection import fit_psp


class TableGroup(object):
    def __init__(self):
        self.mappings = {}
        self.create_mappings()

    def __getitem__(self, item):
        return self.mappings[item]

    def create_mappings(self):
        for k,schema in self.schemas.items():
            self.mappings[k] = db.generate_mapping(k, schema)

    def drop_tables(self):
        for k in self.schemas:
            if k in db.engine.table_names():
                self[k].__table__.drop(bind=db.engine)

    def create_tables(self):
        for k in self.schemas:
            if k not in db.engine.table_names():
                self[k].__table__.create(bind=db.engine)


class PulseResponseStrengthTableGroup(TableGroup):
    """Measures pulse amplitudes for each pulse response and background chunk.
    """
    schemas = {
        'pulse_response_strength': [
            ('pulse_response_id', 'pulse_response.id', '', {'index': True}),
            ('pos_amp', 'float'),
            ('neg_amp', 'float'),
            ('pos_dec_amp', 'float'),
            ('neg_dec_amp', 'float'),
            ('pos_dec_latency', 'float'),
            ('neg_dec_latency', 'float'),
            ('crosstalk', 'float'),
        ],
        'baseline_response_strength' : [
            ('baseline_id', 'baseline.id', '', {'index': True}),
            ('pos_amp', 'float'),
            ('neg_amp', 'float'),
            ('pos_dec_amp', 'float'),
            ('neg_dec_amp', 'float'),
            ('pos_dec_latency', 'float'),
            ('neg_dec_latency', 'float'),
            ('crosstalk', 'float'),
        ]
        #'deconv_pulse_response': [
            #"Exponentially deconvolved pulse responses",
        #],
    }

    def create_mappings(self):
        TableGroup.create_mappings(self)
        
        PulseResponseStrength = self['pulse_response_strength']
        BaselineResponseStrength = self['baseline_response_strength']
        
        db.PulseResponse.pulse_response_strength = db.relationship(PulseResponseStrength, back_populates="pulse_response", cascade="delete", single_parent=True, uselist=False)
        PulseResponseStrength.pulse_response = db.relationship(db.PulseResponse, back_populates="pulse_response_strength", single_parent=True)

        db.Baseline.baseline_response_strength = db.relationship(BaselineResponseStrength, back_populates="baseline", cascade="delete", single_parent=True, uselist=False)
        BaselineResponseStrength.baseline = db.relationship(db.Baseline, back_populates="baseline_response_strength", single_parent=True)


class ConnectionStrengthTableGroup(TableGroup):
    schemas = {
        'connection_strength': [
            ('pair_id', 'pair.id', '', {'index': True}),
            ('synapse_type', 'str', '"ex" or "in"'),

            # current clamp metrics
            ('ic_n_samples', 'int'),
            ('ic_crosstalk_mean', 'float'),
            ('ic_base_crosstalk_mean', 'float'),
            # amplitude,
            ('ic_amp_mean', 'float'),
            ('ic_amp_stdev', 'float'),
            ('ic_base_amp_mean', 'float'),
            ('ic_base_amp_stdev', 'float'),
            ('ic_amp_ttest', 'float'),
            ('ic_amp_ks2samp', 'float'),
            # deconvolved amplitide
            ('ic_deconv_amp_mean', 'float'),
            ('ic_deconv_amp_stdev', 'float'),
            ('ic_base_deconv_amp_mean', 'float'),
            ('ic_base_deconv_amp_stdev', 'float'),
            ('ic_deconv_amp_ttest', 'float'),
            ('ic_deconv_amp_ks2samp', 'float'),
            # latency
            ('ic_latency_mean', 'float'),
            ('ic_latency_stdev', 'float'),
            ('ic_base_latency_mean', 'float'),
            ('ic_base_latency_stdev', 'float'),
            ('ic_latency_ttest', 'float'),
            ('ic_latency_ks2samp', 'float'),
            
            # voltage clamp metrics
            ('vc_n_samples', 'int'),
            ('vc_crosstalk_mean', 'float'),
            ('vc_base_crosstalk_mean', 'float'),
            # amplitude,
            ('vc_amp_mean', 'float'),
            ('vc_amp_stdev', 'float'),
            ('vc_base_amp_mean', 'float'),
            ('vc_base_amp_stdev', 'float'),
            ('vc_amp_ttest', 'float'),
            ('vc_amp_ks2samp', 'float'),
            # deconvolved amplitide
            ('vc_deconv_amp_mean', 'float'),
            ('vc_deconv_amp_stdev', 'float'),
            ('vc_base_deconv_amp_mean', 'float'),
            ('vc_base_deconv_amp_stdev', 'float'),
            ('vc_deconv_amp_ttest', 'float'),
            ('vc_deconv_amp_ks2samp', 'float'),
            # latency
            ('vc_latency_mean', 'float'),
            ('vc_latency_stdev', 'float'),
            ('vc_base_latency_mean', 'float'),
            ('vc_base_latency_stdev', 'float'),
            ('vc_latency_ttest', 'float'),
            ('vc_latency_ks2samp', 'float'),

            # Average pulse responses
            ('ic_average_response', 'array'),
            ('ic_average_response_t0', 'float'),
            ('ic_average_base_stdev', 'float'),
            ('vc_average_response', 'array'),
            ('vc_average_response_t0', 'float'),
            ('vc_average_base_stdev', 'float'),

            # PSP fit parameters
            ('ic_fit_amp', 'float'),
            ('ic_fit_xoffset', 'float'),
            ('ic_fit_yoffset', 'float'),
            ('ic_fit_rise_time', 'float'),
            ('ic_fit_rise_power', 'float'),
            ('ic_fit_decay_tau', 'float'),
            ('ic_fit_exp_amp', 'float'),
            ('ic_fit_nrmse', 'float'),

            ('vc_fit_amp', 'float'),
            ('vc_fit_xoffset', 'float'),
            ('vc_fit_yoffset', 'float'),
            ('vc_fit_rise_time', 'float'),
            ('vc_fit_rise_power', 'float'),
            ('vc_fit_decay_tau', 'float'),
            ('vc_fit_exp_amp', 'float'),
            ('vc_fit_nrmse', 'float'),

        ],
    }

    def create_mappings(self):
        TableGroup.create_mappings(self)
        
        ConnectionStrength = self['connection_strength']
        
        db.Pair.connection_strength = db.relationship(ConnectionStrength, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
        ConnectionStrength.pair = db.relationship(db.Pair, back_populates="connection_strength", single_parent=True)


pulse_response_strength_tables = PulseResponseStrengthTableGroup()
connection_strength_tables = ConnectionStrengthTableGroup()
PulseResponseStrength = pulse_response_strength_tables['pulse_response_strength']
BaselineResponseStrength = pulse_response_strength_tables['baseline_response_strength']
ConnectionStrength = connection_strength_tables['connection_strength']


def init_tables():
    global PulseResponseStrength, BaselineResponseStrength, ConnectionStrength
    pulse_response_strength_tables.create_tables()
    connection_strength_tables.create_tables()

    PulseResponseStrength = pulse_response_strength_tables['pulse_response_strength']
    BaselineResponseStrength = pulse_response_strength_tables['baseline_response_strength']
    ConnectionStrength = connection_strength_tables['connection_strength']


def measure_peak(trace, sign, spike_time, pulse_times, spike_delay=1e-3, response_window=4e-3):
    # Start measuring response after the pulse has finished, and no earlier than 1 ms after spike onset
    # response_start = max(spike_time + spike_delay, pulse_times[1])

    # Start measuring after spike and hope that the pulse offset doesn't get in the way
    # (if we wait for the pulse to end, then we miss too many fast rise / short latency events)
    response_start = spike_time + spike_delay
    response_stop = response_start + response_window

    # measure baseline from beginning of data until 50µs before pulse onset
    baseline_start = 0
    baseline_stop = pulse_times[0] - 50e-6

    baseline = float_mode(trace.time_slice(baseline_start, baseline_stop).data)
    response = trace.time_slice(response_start, response_stop)

    if sign == '+':
        i = np.argmax(response.data)
    else:
        i = np.argmin(response.data)
    peak = response.data[i]
    latency = response.time_values[i] - spike_time
    return peak - baseline, latency


def measure_sum(trace, sign, baseline=(0e-3, 9e-3), response=(12e-3, 17e-3)):
    baseline = trace.time_slice(*baseline).data.sum()
    peak = trace.time_slice(*response).data.sum()
    return peak - baseline

        
def deconv_filter(trace, pulse_times, tau=15e-3, lowpass=1000., lpf=True, remove_artifacts=False, bsub=True):
    dec = exp_deconvolve(trace, tau)

    if remove_artifacts:
        # after deconvolution, the pulse causes two sharp artifacts; these
        # must be removed before LPF
        cleaned = remove_crosstalk_artifacts(dec, pulse_times)
    else:
        cleaned = dec

    if bsub:
        baseline = np.median(cleaned.time_slice(cleaned.t0+5e-3, cleaned.t0+10e-3).data)
        b_subbed = cleaned-baseline
    else:
        b_subbed = cleaned

    if lpf:
        return filter.bessel_filter(b_subbed, lowpass)
    else:
        return b_subbed


def remove_crosstalk_artifacts(data, pulse_times):
    dt = data.dt
    r = [-50e-6, 250e-6]
    edges = [(int((t+r[0])/dt), int((t+r[1])/dt)) for t in pulse_times]
    # If window is too shortm then it becomes seneitive to sample noise.
    # If window is too long, then it becomes sensitive to slower signals (like the AP following pulse onset)
    return filter.remove_artifacts(data, edges, window=100e-6)


@db.default_session
def rebuild_strength(parallel=True, workers=6, session=None):
    for source in ['pulse_response', 'baseline']:
        print("Rebuilding %s strength table.." % source)
        
        # Divide workload into equal parts for each worker
        max_pulse_id = session.execute('select max(id) from %s' % source).fetchone()[0]
        chunk = 1 + (max_pulse_id // workers)
        parts = [(source, chunk*i, chunk*(i+1)) for i in range(workers)]

        if parallel:
            pool = multiprocessing.Pool(processes=workers)
            pool.map(compute_strength, parts)
        else:
            for part in parts:
                compute_strength(part)


def compute_strength(inds, session=None):
    # Thin wrapper just to allow calling from pool.map
    return _compute_strength(inds, session=session)


def response_query(session):
    """
    Build a query to get all pulse responses along with presynaptic pulse and spike timing
    """
    q = session.query(
        db.PulseResponse.id.label('response_id'),
        db.PulseResponse.data,
        db.PulseResponse.start_time.label('rec_start'),
        db.StimPulse.onset_time.label('pulse_start'),
        db.StimPulse.duration.label('pulse_dur'),
        db.StimSpike.max_dvdt_time.label('spike_time'),
        db.PatchClampRecording.clamp_mode,
        db.PulseResponse.ex_qc_pass,
        db.PulseResponse.in_qc_pass,
    )
    q = q.join(db.StimPulse)
    q = q.join(db.StimSpike)
    q = q.join(db.PulseResponse.recording).join(db.PatchClampRecording) 

    # return qc-failed records as well so we can verify qc is working
    #q = q.filter(((db.PulseResponse.ex_qc_pass==True) | (db.PulseResponse.in_qc_pass==True)))

    return q


def baseline_query(session):
    """
    Build a query to get all baseline responses
    """
    q = session.query(
        db.Baseline.id.label('response_id'),
        db.Baseline.data,
        db.PatchClampRecording.clamp_mode,
        db.Baseline.ex_qc_pass,
        db.Baseline.in_qc_pass,
    )
    q = q.join(db.Baseline.recording).join(db.PatchClampRecording) 

    # return qc-failed records as well so we can verify qc is working
    # q = q.filter(((db.Baseline.ex_qc_pass==True) | (db.Baseline.in_qc_pass==True)))

    return q


def analyze_response_strength(rec, source, remove_artifacts=False, lpf=True, bsub=True, lowpass=1000):
    """Perform a standardized strength analysis on a record selected by response_query or baseline_query.

    1. Determine timing of presynaptic stimulus pulse edges and spike
    2. Measure peak deflection on raw trace
    3. Apply deconvolution / artifact removal / lpf
    4. Measure peak deflection on deconvolved trace
    """
    data = Trace(rec.data, sample_rate=db.default_sample_rate)
    if source == 'pulse_response':
        # Find stimulus pulse edges for artifact removal
        start = rec.pulse_start - rec.rec_start
        pulse_times = [start, start + rec.pulse_dur]
        if rec.spike_time is None:
            # these pulses failed QC, but we analyze them anyway to make all data visible
            spike_time = 11e-3
        else:
            spike_time = rec.spike_time - rec.rec_start
    elif source == 'baseline':
        # Fake stimulus information to ensure that background data receives
        # the same filtering / windowing treatment
        pulse_times = [10e-3, 12e-3]
        spike_time = 11e-3
    else:
        raise ValueError("Invalid source %s" % source)

    results = {}

    results['raw_trace'] = data
    results['pulse_times'] = pulse_times
    results['spike_time'] = spike_time

    # Measure crosstalk from pulse onset
    p1 = data.time_slice(pulse_times[0]-200e-6, pulse_times[0]).median()
    p2 = data.time_slice(pulse_times[0], pulse_times[0]+200e-6).median()
    results['crosstalk'] = p2 - p1

    # crosstalk artifacts in VC are removed before deconvolution
    if rec.clamp_mode == 'vc' and remove_artifacts is True:
        data = remove_crosstalk_artifacts(data, pulse_times)
        remove_artifacts = False

    # Measure deflection on raw data
    results['pos_amp'], _ = measure_peak(data, '+', spike_time, pulse_times)
    results['neg_amp'], _ = measure_peak(data, '-', spike_time, pulse_times)

    # Deconvolution / artifact removal / filtering
    tau = 15e-3 if rec.clamp_mode == 'ic' else 5e-3
    dec_data = deconv_filter(data, pulse_times, tau=tau, lpf=lpf, remove_artifacts=remove_artifacts, bsub=bsub, lowpass=lowpass)
    results['dec_trace'] = dec_data

    # Measure deflection on deconvolved data
    results['pos_dec_amp'], results['pos_dec_latency'] = measure_peak(dec_data, '+', spike_time, pulse_times)
    results['neg_dec_amp'], results['neg_dec_latency'] = measure_peak(dec_data, '-', spike_time, pulse_times)
    
    return results


@db.default_session
def _compute_strength(inds, session=None):
    """Comput per-pulse-response strength metrics
    """
    source, start_id, stop_id = inds
    if source == 'baseline':
        q = baseline_query(session)
    elif source == 'pulse_response':
        q = response_query(session)
    else:
        raise ValueError("Invalid source %s" % source)

    prof = pg.debug.Profiler(delayed=False)
    
    next_id = start_id
    while True:
        # Request just a chunk of all pulse responses based on ID range
        if source == 'pulse_response':
            q1 = q.filter(db.PulseResponse.id>=next_id).filter(db.PulseResponse.id<stop_id)
            q1 = q1.order_by(db.PulseResponse.id)
        else:
            q1 = q.filter(db.Baseline.id>=next_id).filter(db.Baseline.id<stop_id)
            q1 = q1.order_by(db.Baseline.id)
        q1 = q1.limit(1000)  # process in 1000-record chunks

        prof('exec')
        recs = q1.all()
        prof('fetch')
        if len(recs) == 0:
            break
        new_recs = []

        for rec in recs:
            result = analyze_response_strength(rec, source)
            new_rec = {'%s_id'%source: rec.response_id}
            # copy a subset of results over to new record
            for k in ['pos_amp', 'neg_amp', 'pos_dec_amp', 'neg_dec_amp', 'pos_dec_latency', 'neg_dec_latency', 'crosstalk']:
                new_rec[k] = result[k]
            new_recs.append(new_rec)
        
        next_id = rec.response_id + 1
        
        sys.stdout.write("%d / %d\r" % (next_id-start_id, stop_id-start_id))
        sys.stdout.flush()

        prof('process')
        if source == 'pulse_response':
            session.bulk_insert_mappings(PulseResponseStrength, new_recs)
        else:
            session.bulk_insert_mappings(BaselineResponseStrength, new_recs)
        prof('insert')
        new_recs = []

        session.commit()
        prof('commit')


@db.default_session
def rebuild_connectivity(session):
    print("Rebuilding connectivity table..")
    
    expts_in_db = list_experiments(session=session)
    for i,expt in enumerate(expts_in_db):
        for pair in expt.pairs:
            # Query all pulse amplitude records for each clamp mode
            amps = {}
            for clamp_mode in ('ic', 'vc'):
                clamp_mode_fg = get_amps(session, pair, clamp_mode=clamp_mode, get_data=True)
                clamp_mode_bg = get_baseline_amps(session, pair, amps=clamp_mode_fg, clamp_mode=clamp_mode, get_data=False)
                amps[clamp_mode, 'fg'] = clamp_mode_fg
                amps[clamp_mode, 'bg'] = clamp_mode_bg
            
            # Generate summary results for this pair
            results = analyze_pair_connectivity(amps)
            
            # Write new record to DB
            conn = ConnectionStrength(pair_id=pair.id, **results)
            session.add(conn)
        
        session.commit()
        sys.stdout.write("%d / %d       \r" % (i, len(expts_in_db)))
        sys.stdout.flush()


def norm_pvalue(pval):
    """Normalize a p-value into a nice 0-7ish range.

    Typically p values can have very large negative exponents, which
    makes them difficult to visualize and to use as classifier features.
    This function is a monotonic remapping of the entire 64-bit floating-point
    space from (~2.5e-324, 1.0) onto a more evenly distributed scale (7, 0).
    """
    return min(7, np.log(1-np.log(pval)))


def analyze_pair_connectivity(amps, sign=None):
    """Given response strength records for a single pair, generate summary
    statistics characterizing strength, latency, and connectivity.
    
    Parameters
    ----------
    amps : dict
        Contains foreground and background strength analysis records
        (see input format below)
    sign : None, -1, or +1
        If None, then automatically determine whether to treat this connection as
        inhibitory or excitatory.

    Input must have the following structure::
    
        amps = {
            ('ic', 'fg'): recs, 
            ('ic', 'bg'): recs,
            ('vc', 'fg'): recs, 
            ('vc', 'bg'): recs,
        }
        
    Where each *recs* must be a structured array containing fields as returned
    by get_amps() and get_baseline_amps().
    
    The overall strategy here is:
    
    1. Make an initial decision on whether to treat this pair as excitatory or
       inhibitory, based on differences between foreground and background amplitude
       measurements
    2. Generate mean and stdev for amplitudes, deconvolved amplitudes, and deconvolved
       latencies
    3. Generate KS test p values describing the differences between foreground
       and background distributions for amplitude, deconvolved amplitude, and
       deconvolved latency    
    """
    requested_sign = sign
    fields = {}  # used to fill the new DB record
    
    # Use KS p value to check for differences between foreground and background
    qc_amps = {}
    ks_pvals = {}
    amp_means = {}
    amp_diffs = {}
    for clamp_mode in ('ic', 'vc'):
        clamp_mode_fg = amps[clamp_mode, 'fg']
        clamp_mode_bg = amps[clamp_mode, 'bg']
        if (len(clamp_mode_fg) == 0 or len(clamp_mode_bg) == 0):
            continue
        for sign in ('pos', 'neg'):
            # Separate into positive/negative tests and filter out responses that failed qc
            qc_field = {'vc': {'pos': 'in_qc_pass', 'neg': 'ex_qc_pass'}, 'ic': {'pos': 'ex_qc_pass', 'neg': 'in_qc_pass'}}[clamp_mode][sign]
            fg = clamp_mode_fg[clamp_mode_fg[qc_field]]
            bg = clamp_mode_bg[clamp_mode_bg[qc_field]]
            qc_amps[sign, clamp_mode, 'fg'] = fg
            qc_amps[sign, clamp_mode, 'bg'] = bg
            if (len(fg) == 0 or len(bg) == 0):
                continue
            
            # Measure some statistics from these records
            fg = fg[sign + '_dec_amp']
            bg = bg[sign + '_dec_amp']
            pval = scipy.stats.ks_2samp(fg, bg).pvalue
            ks_pvals[(sign, clamp_mode)] = pval
            # we could ensure that the average amplitude is in the right direction:
            fg_mean = np.mean(fg)
            bg_mean = np.mean(bg)
            amp_means[sign, clamp_mode] = {'fg': fg_mean, 'bg': bg_mean}
            amp_diffs[sign, clamp_mode] = fg_mean - bg_mean

    if requested_sign is None:
        # Decide whether to treat this connection as excitatory or inhibitory.
        #   strategy: accumulate evidence for either possibility by checking
        #   the ks p-values for each sign/clamp mode and the direction of the deflection
        is_exc = 0
        # print(expt.acq_timestamp, pair.pre_cell.ext_id, pair.post_cell.ext_id)
        for sign in ('pos', 'neg'):
            for mode in ('ic', 'vc'):
                ks = ks_pvals.get((sign, mode), None)
                if ks is None:
                    continue
                # turn p value into a reasonable scale factor
                ks = norm_pvalue(ks)
                dif_sign = 1 if amp_diffs[sign, mode] > 0 else -1
                if mode == 'vc':
                    dif_sign *= -1
                is_exc += dif_sign * ks
                # print("    ", sign, mode, is_exc, dif_sign * ks)
    else:
        is_exc = requested_sign

    if is_exc > 0:
        fields['synapse_type'] = 'ex'
        signs = {'ic':'pos', 'vc':'neg'}
    else:
        fields['synapse_type'] = 'in'
        signs = {'ic':'neg', 'vc':'pos'}

    # compute the rest of statistics for only positive or negative deflections
    for clamp_mode in ('ic', 'vc'):
        sign = signs[clamp_mode]
        fg = qc_amps.get((sign, clamp_mode, 'fg'))
        bg = qc_amps.get((sign, clamp_mode, 'bg'))
        if fg is None or bg is None or len(fg) == 0 or len(bg) == 0:
            fields[clamp_mode + '_n_samples'] = 0
            continue
        
        fields[clamp_mode + '_n_samples'] = len(fg)
        fields[clamp_mode + '_crosstalk_mean'] = np.mean(fg['crosstalk'])
        fields[clamp_mode + '_base_crosstalk_mean'] = np.mean(bg['crosstalk'])
        
        # measure mean, stdev, and statistical differences between
        # fg and bg for each measurement
        for val, field in [('amp', 'amp'), ('deconv_amp', 'dec_amp'), ('latency', 'dec_latency')]:
            f = fg[sign + '_' + field]
            b = bg[sign + '_' + field]
            fields[clamp_mode + '_' + val + '_mean'] = np.mean(f)
            fields[clamp_mode + '_' + val + '_stdev'] = np.std(f)
            fields[clamp_mode + '_base_' + val + '_mean'] = np.mean(b)
            fields[clamp_mode + '_base_' + val + '_stdev'] = np.std(b)
            # statistical tests comparing fg vs bg
            # Note: we use log(1-log(pval)) because it's nicer to plot and easier to
            # use as a classifier input
            tt_pval = scipy.stats.ttest_ind(f, b, equal_var=False).pvalue
            ks_pval = scipy.stats.ks_2samp(f, b).pvalue
            fields[clamp_mode + '_' + val + '_ttest'] = norm_pvalue(tt_pval)
            fields[clamp_mode + '_' + val + '_ks2samp'] = norm_pvalue(ks_pval)


        ### generate the average response and psp fit
        
        # collect all bg and fg traces
        # bg_traces = TraceList([Trace(data, sample_rate=db.default_sample_rate) for data in amps[clamp_mode, 'bg']['data']])
        fg_traces = TraceList()
        for rec in fg:
            t0 = rec['response_start_time'] - rec['max_dvdt_time']   # time-align to presynaptic spike
            trace = Trace(rec['data'], sample_rate=db.default_sample_rate, t0=t0)
            fg_traces.append(trace)
        
        # get averages
        # bg_avg = bg_traces.mean()        
        fg_avg = fg_traces.mean()
        base_rgn = fg_avg.time_slice(-6e-3, 0)
        base = float_mode(base_rgn.data)
        fields[clamp_mode + '_average_response'] = fg_avg.data
        fields[clamp_mode + '_average_response_t0'] = fg_avg.t0
        fields[clamp_mode + '_average_base_stdev'] = base_rgn.std()

        sign = {'pos':'+', 'neg':'-'}[signs[clamp_mode]]
        fg_bsub = fg_avg.copy(data=fg_avg.data - base)  # remove base to help fitting
        try:
            fit = fit_psp(fg_bsub, mode=clamp_mode, sign=sign, xoffset=(1e-3, 0, 6e-3), yoffset=0, mask_stim_artifact=False, rise_time_mult_factor=4)              
            for param, val in fit.best_values.items():
                fields['%s_fit_%s' % (clamp_mode, param)] = val
            fields[clamp_mode + '_fit_yoffset'] = fit.best_values['yoffset'] + base
            fields[clamp_mode + '_fit_nrmse'] = fit.nrmse()
        except:
            print("Error in PSP fit:")
            sys.excepthook(*sys.exc_info())
            continue
        
        #global fit_plot
        #if fit_plot is None:
            #fit_plot = FitExplorer(fit)
            #fit_plot.show()
        #else:
            #fit_plot.set_fit(fit)
        #raw_input("Waiting to continue..")

    return fields
