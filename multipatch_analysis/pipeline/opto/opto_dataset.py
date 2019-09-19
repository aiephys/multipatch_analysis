from collections import OrderedDict
import os
from ... import config #, synphys_cache, lims, qc
from ...util import timestamp_to_datetime
from ..pipeline_module import DatabasePipelineModule
from .opto_experiment import OptoExperimentPipelineModule
from neuroanalysis.data.experiment import Experiment
from neuroanalysis.data.libraries import opto


class OptoDatasetPipelineModule(DatabasePipelineModule):

    name='opto_dataset'
    dependencies = [OptoExperimentPipelineModule]
    table_group = [
        'sync_rec', 
        'recording', 
        'patch_clamp_recording', 
        'test_pulse', 
        'optical_stim_pulse', 
        'pulse_response',
        'baseline'
        ]

    # datasets are large and NWB access leaks memory -- this is probably true for optoanalysis too
    # when running parallel, each child process may run only one job before being killed
    maxtasksperchild = 1 

    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        job_id = job['job_id']

        # Load experiment from DB
        expt_entry = db.experiment_from_ext_id(job_id, session=session)

        # load NWB file
        path = os.path.join(config.synphys_data, expt_entry.storage_path)
        expt = Experiment(site_path=path, loading_library=opto)
        nwb = expt.data

        if expt.ephys_file is not None and expt.uid != '2019_06_13_exp1_TH':
            nwb.contents
            raise Exception('stop')

    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        # only need to return from syncrec table; other tables will be dropped automatically.
        db = self.database
        return session.query(db.SyncRec).filter(db.SyncRec.experiment_id==db.Experiment.id).filter(db.Experiment.ext_id.in_(job_ids)).all()

    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        db = self.database
        
        # All experiments and their creation times in the DB
        expt_module = self.pipeline.get_module('opto_experiment')
        expts = expt_module.finished_jobs()
        
        # Look up nwb file locations for all experiments
        session = db.session()
        expt_recs = session.query(db.Experiment.ext_id, db.Experiment.storage_path, db.Experiment.ephys_file).filter(db.Experiment.ephys_file != None).all()
        expt_paths = {rec.ext_id: rec for rec in expt_recs}
        session.rollback()
        
        # Return the greater of NWB mod time and experiment DB record mtime
        ready = OrderedDict()
        for expt_id, (expt_mtime, expt_success) in expts.items():
            if expt_id not in expt_paths or expt_success is False:
                # no NWB file; ignore
                continue
            rec = expt_paths[expt_id]
            ephys_file = os.path.join(config.synphys_data, rec.storage_path, rec.ephys_file)
            nwb_mtime = timestamp_to_datetime(os.stat(ephys_file).st_mtime)
            ready[rec.ext_id] = max(expt_mtime, nwb_mtime)
        return ready