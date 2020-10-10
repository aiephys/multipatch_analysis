"""
For generating DB tables describing 1) a cells location within cortex, and adding distance info to pairs

"""
from ..pipeline_module import DatabasePipelineModule
from .experiment import ExperimentPipelineModule
from aisynphys import lims
from aisynphys.layer_depths import get_depths_slice
import logging
import numpy as np


class CortexLocationPipelineModule(DatabasePipelineModule):
    """Imports cell location data for each experiment
    """
    name = 'cortical_location'
    dependencies = [ExperimentPipelineModule]
    table_group = ['cortical_cell_location']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        expt_id = job['job_id'] ## an expt id

        expt = db.experiment_from_ext_id(expt_id, session=session)
        slice_entry = expt.slice

        image_series_id = lims.image_series_id(slice_entry.lims_specimen_name)
        # image_series_id = next(image.get('image_series') for image in images if image.get('treatment')=='DAPI')
        cell_specimen_ids = [cell.meta.get('lims_specimen_id') for cell in expt.cell_list]
        cell_specimen_ids = [x for x in cell_specimen_ids if x is not None]
        
        results = get_depths_slice(image_series_id, cell_specimen_ids)
        mean_results = results['mean']
        pia_direction=mean_results['pia_direction']

        for cell in expt.cell_list:
            specimen_id = cell.meta.get('lims_specimen_id')
            if specimen_id not in results:
                continue
            cell_results = results[specimen_id]
            loc_entry = db.CorticalCellLocation(
                layer=cell_results["layer"],
                distance_to_pia=cell_results["absolute_depth"],
                distance_to_wm=cell_results["wm_distance"],
                fractional_depth=cell_results["normalized_depth"],
                layer_depth=cell_results["layer_depth"],
                fractional_layer_depth=cell_results["normalized_layer_depth"],
                position=list(cell_results["position"]),
                cell=cell,
            )
            session.add(loc_entry)

        for pair in expt.pair_list:
            try:
                d12_lat, d12_vert = get_pair_distances(pair, pia_direction)
                pair.lateral_distance = d12_lat
                pair.vertical_distance = d12_vert
            except Exception:
                logging.error("Pair distances not calculated", exc_info=True)

    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        return (session.query(db.CorticalCellLocation)
            .join(db.Cell)
            .join(db.Experiment)
            .filter(db.Experiment.ext_id.in_(job_ids))
            .all())

def get_pair_distances(pair, pia_direction):
    l1 = np.array(pair.pre_cell.cortical_location.position)
    l2 = np.array(pair.post_cell.cortical_location.position)
    if l1 is None or l2 is None:
        raise ValueError(f"Cell locations not found for pair {pair}")
    d12 = l1 - l2
    d12_vert = np.abs(np.dot(d12, pia_direction))[0]
    d12_lat = np.sqrt(np.sum(d12**2) - d12_vert**2)[0]
    return d12_lat, d12_vert