# coding: utf8
"""
For generating a DB table describing short term dynamics.

"""
from __future__ import print_function, division

import os, logging
from collections import OrderedDict
import numpy as np
from ...dynamics import generate_pair_dynamics
from .pipeline_module import MultipatchPipelineModule
from .experiment import ExperimentPipelineModule
from ...stochastic_release_model import list_cached_results, load_cache_file, StochasticReleaseModel
from ...util import datetime_to_timestamp, timestamp_to_datetime
from .dynamics import generate_pair_dynamics


class SynapseModelPipelineModule(MultipatchPipelineModule):
    """Summarizes stochastic model outputs for each synapse
    """
    name = 'synapse_model'
    dependencies = [ExperimentPipelineModule]
    table_group = ['synapse_model']
    
    @classmethod
    def create_db_entries(cls, job, session):
        logger = logging.getLogger(__name__)
        db = job['database']
        job_id = job['job_id']
        if job['meta'] is not None and 'source' in job['meta']:
            cache_file = job['meta']['source']
        else:
            cache_file = cached_results()[job_id][0]
        
        logger.debug("Processing job %s", job_id)

        entry = make_model_result_entry(job_id, db, session, cache_file)

        session.add(entry)
        logger.debug(f"Finished synapse_model for pair {job_id}")
        session.commit()
                    
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        all_recs = session.query(db.SynapseModel).all()
        selected_recs = [r for r in all_recs if r.meta['pair_ext_id'] in job_ids]         
        return selected_recs

    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        logger = logging.getLogger(__name__)

        # find all cached model results and check mtime
        results = cached_results()
        logger.info(f"Found {len(results)} model cache files")

        # when did each experiment pipeline job finish? (when pair entries were created)
        db = self.database
        pair_pipeline_jobs = db.query(db.Pipeline).filter(db.Pipeline.module_name=='experiment').all()
        pair_job_timestamps = {job.job_id: job.finish_time for job in pair_pipeline_jobs if job.success}

        # get pair IDs of all known synapses
        q = db.pair_query(synapse=True)
        q = q.add_columns(
            db.Experiment.ext_id.label('expt_ext_id'),
            q.pre_cell.ext_id.label('pre_ext_id'),
            q.post_cell.ext_id.label('post_ext_id'),
        )
        synapses = q.all()
        pair_ids = set([" ".join([syn.expt_ext_id, syn.pre_ext_id, syn.post_ext_id]) for syn in synapses])

        ready = OrderedDict()
        for pair_id, (cache_file, mtime) in results.items():
            if pair_id not in pair_ids:
                logger.warn(f"Skipping cached model result for pair {pair_id} not in DB ({cache_file})")
                continue

            # take the latter of cache file mtime and the pipeline run time for the pair
            expt_id, _, _ = pair_id.split(' ')
            dep_time = max(timestamp_to_datetime(mtime), pair_job_timestamps[expt_id] ) 

            ready[pair_id] = {'dep_time': dep_time, 'meta': {'source': cache_file}}
        return ready


_all_results = None
def cached_results():
    """Return a dict mapping {pair_id: (cache_file, mtime)} for all cached model results.
    """
    global _all_results
    if _all_results is not None:
        return _all_results
        
    _all_results = {}
    for pair_id, cache_file in list_cached_results():
        pair_id = ' '.join(pair_id)
        mtime = os.stat(cache_file).st_mtime
        _all_results[pair_id] = (cache_file, mtime)
    
    return _all_results


def make_model_result_entry(pair_id, db, session, model_cache_file):
    expt_id, pre_id, post_id = pair_id.split(' ')

    # Load experiment from DB
    expt = db.experiment_from_ext_id(expt_id, session=session)
    pair = expt.pairs[pre_id, post_id]

    # load cached model result
    model_runner = load_cache_file(model_cache_file, db)
    (spikes, amps, baseline, extra) = model_runner.synapse_events
    n_events = int((np.isfinite(amps) & np.isfinite(spikes)).sum())

    entry = db.SynapseModel(
        pair_id=pair.id,
        n_source_events=n_events,
        meta={
            'pair_ext_id': pair_id,
        }
    )


    # get likelihood and amplitude values across the parameter space,
    param_space = model_runner.param_space
    entry.parameter_space = [(k, list(v)) for k,v in param_space.axes().items()]

    max_likelihood, ml_index, ml_params = get_max_likelihood(param_space)
    entry.max_likelihood = max_likelihood
    for k,v in ml_params.items():
        setattr(entry, 'ml_'+k, v)

    strength, quanta_per_spike, sites_pr_ratio = get_model_metrics(ml_index, model_runner)
    entry.ml_strength = strength
    entry.ml_quanta_per_spike = quanta_per_spike
    entry.ml_sites_pr_ratio = sites_pr_ratio

    entry.ml_release_dependence_ratio = get_release_dependence_ratio(param_space)

    # Simulate an experiment and pass the data through the dynamics pipeline to generate similar results
    dynamics = get_model_dynamics(ml_index, model_runner, pair, db)
    # transfer results over to the new entry
    for name in dir(dynamics):
        if hasattr(entry, 'ml_'+name):
            setattr(entry, 'ml_'+name, getattr(dynamics, name))

    return entry


def get_release_dependence_ratio(param_space):
    """Return the ratio of max likelihoods for the release-dependent and
    release-independent portions of the parameter space
    """
    # which result axis is depression_amount?
    dep_axis = list(param_space.axes().keys()).index('depression_amount')

    # separate the release-dependent results from release-independent
    rel_dep_slice = [slice(None)] * param_space.result.ndim
    rel_indep_slice = [slice(None)] * param_space.result.ndim
    rel_dep_slice[dep_axis] = slice(0, 1)
    rel_indep_slice[dep_axis] = slice(1, None)

    rel_dep_result = param_space.result[tuple(rel_dep_slice)]
    rel_indep_result = param_space.result[tuple(rel_indep_slice)]

    rel_dep_ml = rel_dep_result['likelihood'].max()
    rel_indep_ml = rel_indep_result['likelihood'].max()

    return rel_dep_ml / rel_indep_ml


def get_max_likelihood(param_space):
    """Return maximum likelihood, parameter space index, and a dictionary of ML parameter values
    """
    res = param_space.result
    max_likelihood = res['likelihood'].max()
    max_index = list(np.unravel_index(np.argmax(res['likelihood']), res['likelihood'].shape))
    max_params = get_params_from_index(param_space, max_index)

    return max_likelihood, max_index, max_params


def get_params_from_index(param_space, index):
    params = {'mini_amplitude': param_space.result['mini_amplitude'][tuple(index)]}
    params.update(param_space.params_at_index(index))
    return params


def get_model_metrics(index, model_runner):
    """Return metrics derived from model parameters
    """
    params = get_params_from_index(model_runner.param_space, index)
    quanta_per_spike = params['base_release_probability'] * params['n_release_sites']
    strength = quanta_per_spike * params['mini_amplitude'] 
    sites_pr_ratio = params['n_release_sites'] / params['base_release_probability']

    return strength, quanta_per_spike, sites_pr_ratio


def get_model_dynamics(index, model_runner, pair, db):
    """Return dynamics metrics generated from simulated events.
    """
    session = db.session()
    model_pair = db.Pair(
        has_synapse=True,
    )
    model_pair.synapse = db.Synapse(
        synapse_type=pair.synapse.synapse_type,
    )
    model_pair.pre_cell = db.Cell(
        ext_id=pair.pre_cell.ext_id,
    )
    model_pair.post_cell = db.Cell(
        ext_id=pair.post_cell.ext_id,
    )
    model_pair.experiment = db.Experiment(
        ext_id='MOCK_'+pair.experiment.ext_id,
    )

    # spike timing for a 12-pulse train
    ind_freq = 50
    rec_delay = 250e-3
    templ = np.arange(12) / 50
    templ[8:] += rec_delay - 1/ind_freq

    # repeat 500x
    n_rep = 500
    iti = 15
    n_pulses = len(templ)
    spike_times = np.empty(n_pulses * n_rep)
    for i in range(n_rep):
        spike_times[i*n_pulses:(i+1)*n_pulses] = templ + i * iti

    pulse_n = (np.arange(n_pulses * n_rep) % n_pulses)  + 1

    # run simulation
    params = get_params_from_index(model_runner.param_space, index)
    model = StochasticReleaseModel(params)
    result = model.run_model(spike_times, amplitudes='random')
    amps = result.result['amplitude']

    class PrRec:
        """Mock result from aisynphys.dynamics.pulse_response_query
        fields: db.PulseResponse, db.PulseResponseFit, db.StimPulse, db.Recording, db.PatchClampRecording, db.MultiPatchProbe, db.Synapse
        """
        def __init__(self, **kwds):
            for k,v in kwds.items():
                setattr(self, k, v)

    recordings = []
    for i in range(n_rep):
        rec = db.Recording(stim_name="mock_mp_probe")
        rec.electrode = db.Electrode(ext_id=1)
        rec.sync_rec = db.SyncRec(ext_id=i)
        rec.sync_rec.experiment = model_pair.experiment
        rec.patch_clamp_recording = db.PatchClampRecording(clamp_mode='ic')
        rec.patch_clamp_recording.multi_patch_probe = db.MultiPatchProbe(induction_frequency=ind_freq, recovery_delay=rec_delay)
        recordings.append(rec)

    pr_recs = []
    for i in range(len(spike_times)):
        pr = db.PulseResponse(
            ex_qc_pass=True,
            in_qc_pass=True,
        )
        pr.pulse_response_fit = db.PulseResponseFit(
            dec_fit_reconv_amp=amps[i],
            baseline_dec_fit_reconv_amp=0,
        )
        pr.stim_pulse = db.StimPulse(
            pulse_number=pulse_n[i],
            previous_pulse_dt=np.inf if i==0 else spike_times[i]-spike_times[i-1],
        )
        pr.recording = recordings[i//n_pulses]

        pr_recs.append(PrRec(
            PulseResponse=pr,
            PulseResponseFit=pr.pulse_response_fit, 
            StimPulse=pr.stim_pulse,
            Recording=pr.recording,
            PatchClampRecording=pr.recording.patch_clamp_recording,
            MultiPatchProbe=pr.recording.patch_clamp_recording.multi_patch_probe,
            Synapse=model_pair.synapse
        ))


    dynamics = generate_pair_dynamics(model_pair, db, session, pr_recs)
    return dynamics
