from aisynphys.pipeline.pipeline_module import DatabasePipelineModule
from .opto_experiment import OptoExperimentPipelineModule
import aisynphys.database as db

class OptoMorphologyPipelineModule(DatabasePipelineModule):
    """Imports cell morphology data for each experiment
    """
    name = 'opto_morphology'
    dependencies = [OptoExperimentPipelineModule]
    #table_group = db.morphology_tables
    table_group = ['morphology']

    @classmethod
    def create_db_entries(cls, job_id, session, expt=None):

        if expt == None:
            raise NotImplementedError("Adding morphology without an Experiment data model object is not yet implemented.")
            expt = db.experiment_from_uid(job_id, session=session)

        db_expt = db.experiment_from_uid(job_id, session=session)

        try:

            for name, cell in expt.cells.items():
                db_cell = db_expt.cells[name]
                morph = cell.morphology.get('initial_call')
                if morph == 'pyramidal':
                    pyr = True
                elif morph == 'interneuron':
                    pyr = False
                else:
                    pyr = None
                #print('cell_entry.id:', cell_entry.id)
                #print('cell.cell_id:', cell.cell_id)
                morphology_entry = db.Morphology(
                    cell_id=db_cell.id,
                    pyramidal=pyr)
                session.add(morphology_entry)

        except:
            session.rollback()

        session.commit()

    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        return session.query(db.Morphology).filter(db.Morphology.cell_id==db.Cell.id).filter(db.Cell.experiment_id==db.Experiment.id).filter(db.Experiment.ext_id.in_(job_ids)).all()

    # def ready_jobs(self):
    #     """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
    #     and the dates that dependencies were created.
    #     """
    #     db = self.database
    #     # All experiments and their creation times in the DB
    #     expts = self.pipeline.get_module('opto_experiment').finished_jobs()

    #     # Look up nwb file locations for all experiments
    #     session = db.session()
    #     # expts = session.query(db.Experiment).filter(db.Experiment.acq_timestamp==1521667891.153).all()
    #     session.rollback()
        
    #     # Return the greater of NWB mod time and experiment DB record mtime
    #     ready = OrderedDict()

    #     try:
    #         morpho_results = morpho_db()
    #     except ImportError as exc:
    #         print("Skipping morphology: %s" % str(exc))
    #         return ready

    #     for expt_id, (expt_mtime, success) in expts.items():
    #         if success is not True:
    #             continue

    #         expt = session.query(db.Experiment).filter(db.Experiment.acq_timestamp==expt_id).all()[0]
    #         cluster = expt.lims_specimen_id
    #         if cluster is None:
    #             continue
    #         cluster_cells = lims.cluster_cells(cluster)
    #         if cluster_cells is None:
    #             continue
    #         cell_hash_compare = []
    #         for cell in cluster_cells:
    #             if cell.id is None:
    #                 continue
    #             morpho_db_hash = hash(';'.join(filter(None, morpho_results.get(cell.id, [None]))))
    #             cell_morpho_rec = session.query(db.Morphology).join(db.Cell).filter(db.Cell.meta.info.get('lims_specimen_id')==str(cell.id)).all()
    #             if len(cell_morpho_rec) == 1:   
    #                 cell_rec_hash = cell_morpho_rec[0].morpho_db_hash
    #                 cell_hash_compare.append(morpho_db_hash == cell_rec_hash)
    #             else:
    #                 cell_hash_compare.append(False)
    #         if all(cell_hash_compare) is False:
    #             ready[expt.acq_timestamp] = datetime.datetime.now() #### if we use this we'll need to update it to {'dep_time': time, 'meta':{}}
    
    #     return ready