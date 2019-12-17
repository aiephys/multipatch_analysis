"""
For generating DB tables describing 1) a cells location within cortex, 2) a cortical site.

"""
from __future__ import print_function, division
from ..pipeline_module import DatabasePipelineModule
from .opto_experiment import OptoExperimentPipelineModule, load_experiment
from collections import OrderedDict
import datetime, os


class OptoCortexLocationPipelineModule(DatabasePipelineModule):
    """Imports cell morphology data for each experiment
    """
    name = 'cortical_location'
    dependencies = [OptoExperimentPipelineModule]
    table_group = ['cell_location', 'cortical_site']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        job_id = job['job_id'] ## an expt id

        expt = load_experiment(job_id)

        try:
            # look up slice record in DB
            try:
                ts = expt.slice_timestamp
            except KeyError:
                ts = 0.0
            slice_entry = db.slice_from_timestamp(ts, session=session)

            site_entry = db.CorticalSite(
                    pia_to_wm_distance=expt.cortical_site_info.get('pia_to_wm_distance'),
                    pia_position=expt.cortical_site_info.get('piaPos'),
                    wm_position=expt.cortical_site_info.get('wmPos'),
                    L1_L23_boundary=expt.cortical_site_info.get('layerBounds_percentDepth',{}).get('L2/3', [None])[0],
                    L23_L4_boundary=expt.cortical_site_info.get('layerBounds_percentDepth',{}).get('L4', [None])[0],
                    L4_L5_boundary=expt.cortical_site_info.get('layerBounds_percentDepth',{}).get('L5', [None])[0],
                    L5_L6_boundary=expt.cortical_site_info.get('layerBounds_percentDepth',{}).get('L6', [None])[0]
                    )

            site_entry.slice = slice_entry
            site_entry.experiment = db.experiment_from_ext_id(expt.uid, session=session)
            session.add(site_entry)


            for cell_id, cell in expt.cells.items():

                loc_entry = db.CellLocation(
                    #cell_id=cell.cell_id,
                    layer=cell.target_layer,
                    distance_to_pia=cell.distance_to_pia,
                    distance_to_wm=cell.distance_to_wm,
                    fractional_depth=cell.percent_depth
                )
                cell_entry = session.query(db.Cell).filter(db.Cell.experiment_id==db.Experiment.id).filter(db.Experiment.ext_id==job_id).filter(db.Cell.ext_id==cell_id).all()
                if len(cell_entry) == 1:
                    loc_entry.cell = cell_entry[0]
                else:
                    raise Exception("Found wrong number of cell entries for experiment %s, cell %s" %(job_id, cell_id))
                loc_entry.site = site_entry
                session.add(loc_entry)

        except:
            session.rollback()
            raise


    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        return session.query(db.CorticalSite).filter(db.CorticalSite.experiment_id==db.Experiment.id).filter(db.Experiment.ext_id.in_(job_ids)).all()

        #return session.query(db.Morphology).filter(db.Morphology.cell_id==db.Cell.id).filter(db.Cell.experiment_id==db.Experiment.id).filter(db.Experiment.acq_timestamp.in_(job_ids)).all()

    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        db = self.database
        # All experiments and their creation times in the DB
        expts = self.pipeline.get_module('opto_experiment').finished_jobs()

        session = db.session()
        session.rollback()
        
        ready = OrderedDict()

        for expt_id, (expt_mtime, success) in expts.items():
            if success is not True:
                continue
            expt = load_experiment(expt_id)
            if expt.connections_file_version <= 3:
                mtime = datetime.datetime.fromtimestamp(os.path.getmtime(expt.connections_file))
            else:
                mtime = datetime.datetime.fromtimestamp(os.path.getmtime(config.distance_csv))
            ready[expt_id] = {'dep_time':mtime, 'meta':{}}
    
        return ready