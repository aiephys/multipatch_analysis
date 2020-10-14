"""
For generating DB table describing a cells location within cortex, 
and adding layer-aligned distance info to the Pair table
"""
from ..pipeline_module import DatabasePipelineModule
from .experiment import ExperimentPipelineModule
from aisynphys import lims
import logging
import numpy as np
from neuroanalysis.util.optional_import import optional_import
get_depths_slice = optional_import('aisynphys.layer_depths', 'get_depths_slice')


class CortexLocationPipelineModule(DatabasePipelineModule):
    """Imports cell location data for each experiment
    """
    name = 'cortical_location'
    dependencies = [ExperimentPipelineModule]
    table_group = ['cortical_cell_location']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        expt_id = job['job_id']

        expt = db.experiment_from_ext_id(expt_id, session=session)
        slice_entry = expt.slice

        try:
            images = lims.specimen_images(slice_entry.lims_specimen_name)
            images = [image for image in images if image.get('treatment')=='DAPI']
            assert len(images) > 0
            image_series_id = images[0].get('image_series')
            image_series_resolution = images[0].get('resolution')
        except (AssertionError, ValueError):
            logging.warning("No LIMS image series found for slice.")
            return

        try:
            lims_cell_cluster_id = expt.meta.get('lims_cell_cluster_id')
            lims_cell_info = lims.cluster_cells(lims_cell_cluster_id)
            soma_centers = {cell['id']: (cell['x_coord'], cell['y_coord']) for cell in lims_cell_info}
            soma_centers = {cell: coords for cell, coords in soma_centers.items() 
                            if all(coords)}
            assert len(soma_centers) > 0
        except (AssertionError, ValueError):
            logging.warning("No cell coordinates found for cell cluster.")
            return
        
        results = get_depths_slice(image_series_id, soma_centers, resolution=image_series_resolution)
        if len(results)==0:
            logging.error("No cells passed depth calculation.")
            return

        missed_cell_count = 0
        for cell in expt.cell_list:
            specimen_id = cell.meta.get('lims_specimen_id')
            if specimen_id not in results:
                missed_cell_count += 1
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
        
        if missed_cell_count > 0:
            logging.warning(f"{missed_cell_count}/{len(expt.cell_list)} cells failed depth calculation.")

        if not 'mean' in results:
            logging.error("Depth calculation for cell cluster mean failed.")
            return
        mean_results = results['mean']
        pia_direction = mean_results['pia_direction']

        for pair in expt.pair_list:
            pre_id = pair.pre_cell.meta.get('lims_specimen_id')
            post_id = pair.post_cell.meta.get('lims_specimen_id')
            if pre_id in results and post_id in results:
                try:
                    d12_lat, d12_vert = get_pair_distances(pair, pia_direction)
                    pair.lateral_distance = d12_lat
                    pair.vertical_distance = d12_vert
                except Exception:
                    logging.error("Pair distances calculation failed", exc_info=True)

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
    d12 = l1 - l2
    d12_vert = np.abs(np.dot(d12, pia_direction))[0]
    d12_lat = np.sqrt(np.sum(d12**2) - d12_vert**2)
    return d12_lat, d12_vert