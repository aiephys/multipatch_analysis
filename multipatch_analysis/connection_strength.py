# coding: utf8
"""
Analyses that measure the strength of synaptic connections.

"""
from __future__ import print_function, division

import sys, multiprocessing

import numpy as np
import pandas
import scipy.stats

from neuroanalysis.data import Trace, TraceList
from neuroanalysis.baseline import float_mode

from .database import database as db
from .pulse_response_strength import PulseResponseStrength, BaselineResponseStrength
from .connection_detection import fit_psp
from .database import TableGroup


class ConnectionStrengthTableGroup(TableGroup):
    schemas = {
        'connection_strength': [
            """Describes the statistics of per-pair properties aggregated from the pulse_response_strength table.
            """,
            ('pair_id', 'pair.id', 'The ID of the entry in the pair table to which these results apply', {'index': True}),
            ('synapse_type', 'str', 'String "ex" or "in", indicating whether this analysis chose to treat the pair as excitatory or inhibitory'),

            # current clamp metrics
            ('ic_n_samples', 'int', "Number of samples (pulse responses) that were pooled from current clamp recordings"),
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


connection_strength_tables = ConnectionStrengthTableGroup()


def init_tables():
    global ConnectionStrength
    connection_strength_tables.create_tables()

    ConnectionStrength = connection_strength_tables['connection_strength']


# create tables in database and add global variables for ORM classes
init_tables()



@db.default_session
def update_connection_strength(limit=0, expts=None, parallel=True, workers=6, raise_exceptions=False, session=None):
    """Update connection strength table for all experiments
    """
    if expts is None:
        expts_ready = session.query(db.Experiment.acq_timestamp).join(db.SyncRec).join(db.Recording).join(db.PulseResponse).join(PulseResponseStrength).distinct().all()
        expts_done = session.query(db.Experiment.acq_timestamp).join(db.Pair).join(ConnectionStrength).distinct().all()

        print("Skipping %d already complete experiments" % (len(expts_done)))
        experiments = [e for e in expts_ready if e not in set(expts_done)]

        if limit > 0:
            np.random.shuffle(experiments)
            experiments = experiments[:limit]

        jobs = [(record.acq_timestamp, index, len(experiments)) for index, record in enumerate(experiments)]
    else:
        jobs = [(expt, i, len(expts)) for i, expt in enumerate(expts)]

    if parallel:
        pool = multiprocessing.Pool(processes=workers)
        pool.map(compute_connection_strength, jobs)
    else:
        for job in jobs:
            compute_connection_strength(job, raise_exceptions=raise_exceptions)


def compute_connection_strength(job_info, raise_exceptions=False):
    session = db.Session()
    
    try:
        expt_id, index, n_jobs = job_info
        print("Analyzing connection strength (expt_id=%f): %d/%d" % (expt_id, index, n_jobs))

        expt = db.experiment_from_timestamp(expt_id, session=session)

        for pair in expt.pairs:
            # Query all pulse amplitude records for each clamp mode
            amps = {}
            for clamp_mode in ('ic', 'vc'):
                clamp_mode_fg = get_amps(session, pair, clamp_mode=clamp_mode, get_data=True)
                clamp_mode_bg = get_baseline_amps(session, pair, amps=clamp_mode_fg, clamp_mode=clamp_mode, get_data=False)
                amps[clamp_mode, 'fg'] = clamp_mode_fg
                amps[clamp_mode, 'bg'] = clamp_mode_bg
            
            if all([len(a) == 0 for a in amps]):
                # nothing to analyze here.
                continue

            # Generate summary results for this pair
            results = analyze_pair_connectivity(amps)

            # Write new record to DB
            conn = ConnectionStrength(pair_id=pair.id, **results)
            session.add(conn)
        session.commit()
    except:
        session.rollback()
        print("Error in experiment: %f" % expt_id)
        if raise_exceptions:
            raise
        else:
            sys.excepthook(*sys.exc_info())


def get_amps(session, pair, clamp_mode='ic', get_data=False):
    """Select records from pulse_response_strength table
    """
    cols = [
        PulseResponseStrength.id,
        PulseResponseStrength.pos_amp,
        PulseResponseStrength.neg_amp,
        PulseResponseStrength.pos_dec_amp,
        PulseResponseStrength.neg_dec_amp,
        PulseResponseStrength.pos_dec_latency,
        PulseResponseStrength.neg_dec_latency,
        PulseResponseStrength.crosstalk,
        db.PulseResponse.ex_qc_pass,
        db.PulseResponse.in_qc_pass,
        db.PatchClampRecording.clamp_mode,
        db.StimPulse.pulse_number,
        db.StimSpike.max_dvdt_time,
        db.PulseResponse.start_time.label('response_start_time'),
    ]
    if get_data:
        cols.append(db.PulseResponse.data)

    q = session.query(*cols)
    q = q.join(db.PulseResponse)
    
    q, pre_rec, post_rec = join_pulse_response_to_expt(q)
    q = q.join(db.StimSpike)
    q = q.add_columns(post_rec.start_time.label('rec_start_time'))
        
    filters = [
        (pre_rec.electrode==pair.pre_cell.electrode,),
        (post_rec.electrode==pair.post_cell.electrode,),
        (db.PatchClampRecording.clamp_mode==clamp_mode,),
        (db.PatchClampRecording.qc_pass==True,),
    ]
    for filter_args in filters:
        q = q.filter(*filter_args)
    
    # should result in chronological order
    q = q.order_by(db.PulseResponse.id)

    df = pandas.read_sql_query(q.statement, q.session.bind)
    recs = df.to_records()
    return recs


def get_baseline_amps(session, pair, clamp_mode='ic', amps=None, get_data=True):
    """Select records from baseline_response_strength table

    If *amps* is given (output from get_amps), then baseline records will be selected from the same
    sweeps as the responses.
    """
    cols = [
        BaselineResponseStrength.id,
        BaselineResponseStrength.pos_amp,
        BaselineResponseStrength.neg_amp,
        BaselineResponseStrength.pos_dec_amp,
        BaselineResponseStrength.neg_dec_amp,
        BaselineResponseStrength.pos_dec_latency,
        BaselineResponseStrength.neg_dec_latency,
        BaselineResponseStrength.crosstalk,
        db.Baseline.ex_qc_pass,
        db.Baseline.in_qc_pass,
        db.PatchClampRecording.clamp_mode,
        db.Recording.start_time.label('rec_start_time'),
        db.Baseline.start_time.label('response_start_time'),
    ]
    if get_data:
        cols.append(db.Baseline.data)
        
    q = session.query(*cols)
    q = q.join(db.Baseline).join(db.Recording).join(db.PatchClampRecording).join(db.SyncRec).join(db.Experiment)
    
    filters = [
        (db.Recording.electrode==pair.post_cell.electrode,),
        (db.PatchClampRecording.clamp_mode==clamp_mode,),
        (db.PatchClampRecording.qc_pass==True,),
    ]
    for filter_args in filters:
        q = q.filter(*filter_args)
    
    # should result in chronological order
    q = q.order_by(db.Recording.start_time)

    # if amps is not None:
    #     q = q.limit(len(amps))

    df = pandas.read_sql_query(q.statement, q.session.bind)
    recs = df.to_records()

    if amps is not None:
        # for each record returned from get_amps, return the nearest baseline record
        mask = np.zeros(len(recs), dtype=bool)
        amp_times = amps['rec_start_time'].astype(float)*1e-9 + amps['response_start_time']
        base_times = recs['rec_start_time'].astype(float)*1e-9 + recs['response_start_time']
        for i in range(len(amps)):
            order = np.argsort(np.abs(base_times - amp_times[i]))
            for j in order:
                if mask[j]:
                    continue
                mask[j] = True
                break
        recs = recs[mask]

    return recs


def join_pulse_response_to_expt(query):
    pre_rec = db.aliased(db.Recording)
    post_rec = db.aliased(db.Recording)
    joins = [
        (post_rec, db.PulseResponse.recording),
        (db.PatchClampRecording,),
        (db.SyncRec,),
        (db.Experiment,),
        (db.StimPulse, db.PulseResponse.stim_pulse),
        (pre_rec, db.StimPulse.recording),
    ]
    for join_args in joins:
        query = query.join(*join_args)

    return query, pre_rec, post_rec



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
            fit = fit_psp(fg_bsub, mode=clamp_mode, sign=sign, xoffset=(1e-3, 0, 6e-3), yoffset=(0, None, None), rise_time_mult_factor=4)              
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
