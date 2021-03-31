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
from ...stochastic_release_model import list_cached_results, load_cache_file
from ...util import datetime_to_timestamp, timestamp_to_datetime


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
        logger.debug("Processing job %s", job_id)

        expt_id, pre_id, post_id = job_id.split(' ')

        # Load experiment from DB
        expt = db.experiment_from_ext_id(expt_id, session=session)
        pair = expt.pairs[pre_id, post_id]

        model_runner = load_cache_file(job['meta']['source'], db)
        (spikes, amps, baseline, extra) = model_runner.synapse_events
        n_events = int((np.isfinite(amps) & np.isfinite(spikes)).sum())

        entry = db.SynapseModel(
            pair_id=pair.id,
            n_source_events=n_events,
            meta={
                'pair_ext_id': job_id,
            }
        )
        session.add(entry)
        logger.debug("Finished synapse_model for pair %s", pair)
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
