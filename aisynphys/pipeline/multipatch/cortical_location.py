"""
For generating DB table describing a cells location within cortex, 
and adding layer-aligned distance info to the Pair table
"""
from ..pipeline_module import DatabasePipelineModule
from .experiment import ExperimentPipelineModule
from aisynphys import lims
import numpy as np
from neuroanalysis.util.optional_import import optional_import
get_depths_slice = optional_import('aisynphys.layer_depths', 'get_depths_slice')
import logging
logger = logging.getLogger(__name__)


class CortexLocationPipelineModule(DatabasePipelineModule):
    """Imports cell location data for each experiment
    """
    name = 'cortical_location'
    dependencies = [ExperimentPipelineModule]
    table_group = ['cortical_cell_location']
    
    @classmethod
    def create_db_entries(cls, job, session):
        lims_layers = get_lims_layers() 
        db = job['database']
        expt_id = job['job_id']
        errors = []

        expt = db.experiment_from_ext_id(expt_id, session=session)
        slice_entry = expt.slice

        try:
            images = lims.specimen_images(slice_entry.lims_specimen_name)
            images = [image for image in images if image.get('treatment')=='DAPI']
            assert len(images) > 0
            image_series_id = images[0].get('image_series')
            image_series_resolution = images[0].get('resolution')
        except (AssertionError, ValueError):
            errors.append("No LIMS image series found for slice.")
            return errors

        try:
            lims_cell_cluster_id = expt.meta.get('lims_cell_cluster_id')
            lims_cell_info = lims.cluster_cells(lims_cell_cluster_id)
            lims_ids = [cell.meta.get('lims_specimen_id') for cell in expt.cell_list]
            soma_centers = {cell['id']: (cell['x_coord'], cell['y_coord']) 
                            for cell in lims_cell_info}
            soma_centers = {cell: coords for cell, coords in soma_centers.items() 
                            if all(coords) and cell in lims_ids}
            assert len(soma_centers) > 0
        except (AssertionError, ValueError):
            errors.append("No cell coordinates found for cell cluster.")
            return errors
        
        results, cell_errors = get_depths_slice(image_series_id, soma_centers,
                                                species=slice_entry.species,
                                                resolution=image_series_resolution)
        if len(results)==0:
            errors.append("No cells passed depth calculation.")
            errors.extend(cell_errors)
            return errors

        missed_cells = []
        for cell in expt.cell_list:
            specimen_id = cell.meta.get('lims_specimen_id')
            lims_layer = lims_layers.get(specimen_id, None)
            meta = {'lims_layer': lims_layer}
            if specimen_id not in soma_centers:
                continue
            if specimen_id not in results:
                missed_cells.append(cell.ext_id)
                loc_entry = db.CorticalCellLocation(
                    cell=cell,
                    position=soma_centers[specimen_id],
                    meta=meta,
                )
            else:
                cell_results = results[specimen_id]
                loc_entry = db.CorticalCellLocation(
                    cortical_layer=cell_results["layer"].replace("Layer",''),
                    distance_to_pia=cell_results.get("absolute_depth", np.nan)*1e-6,
                    distance_to_wm=cell_results.get("wm_distance", np.nan)*1e-6,
                    fractional_depth=cell_results.get("normalized_depth", np.nan),
                    layer_depth=cell_results.get("layer_depth", np.nan)*1e-6,
                    layer_thickness=cell_results.get("layer_thickness", np.nan)*1e-6,
                    fractional_layer_depth=cell_results.get("normalized_layer_depth", np.nan),
                    position=list(cell_results["position"]*1e-6),
                    cell=cell,
                    meta=meta,
                )
            session.add(loc_entry)

        if len(missed_cells) > 0:
            msg = f"Cells {missed_cells} (of {len(soma_centers)}) cells failed depth calculation."
            logger.error(msg)
            errors.append(msg)
            errors.extend(cell_errors)
            
        for pair in expt.pair_list:
            pre_id = pair.pre_cell.meta.get('lims_specimen_id')
            post_id = pair.post_cell.meta.get('lims_specimen_id')
            if pre_id in soma_centers and post_id in soma_centers:
                nan_dir = np.nan*np.ones((2,1))
                pre_dir = results[pre_id]['pia_direction'] if pre_id in results else nan_dir
                post_dir = results[post_id]['pia_direction'] if post_id in results else nan_dir
                pia_direction = np.nanmean(np.stack([pre_dir, post_dir]), axis=0)
                if all(pia_direction == pia_direction): # no nans
                    d12_lat, d12_vert = get_pair_distances(pair, pia_direction)
                    pair.lateral_distance = d12_lat
                    pair.vertical_distance = d12_vert
                
        return errors

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

lims_cache = None
def get_lims_layers():
    global lims_cache
    if lims_cache is None:
        lims_q = lims.all_cell_layers()
        lims_cache = {spec_id:layer.lstrip('Layer') for layer, spec_id in lims_q if layer is not None}

    return lims_cache