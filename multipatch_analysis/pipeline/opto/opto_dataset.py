from ..pipeline_module import DatabasePipelineModule
from .opto_experiment import OptoExperimentPipelineModule
from neuroanalysis.data.Experiment import Experiment
from neuroanalysis.data.libraries import opto


class OptoDatasetPipelineModule(DatabasePipelineModule):

    name='opto_dataset'
    depencencies = [OptoExperimentPipelineModule]
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