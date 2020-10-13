import numpy as np
from neuron_morphology.lims_apical_queries import get_data
from neuron_morphology.transforms.pia_wm_streamlines.calculate_pia_wm_streamlines import get_depth_gradient_field
from neuron_morphology.snap_polygons.__main__ import run_snap_polygons, Parser
from neuron_morphology.snap_polygons.types import ensure_path, ensure_linestring, ensure_polygon
import neuron_morphology.layered_point_depths.__main__ as ld 
from shapely.geometry import Polygon, Point, LineString, LinearRing
import logging
logger = logging.getLogger(__name__)

class LayerDepthError(Exception):
    pass

def get_cell_soma_data(cell_specimen_ids):
    # based on query in lims_apical_queries but remove requirement of reconstruction
    # there is probably a better way to do this, and should be in neuron_morphology
    ids_str = ', '.join([str(sid) for sid in cell_specimen_ids])
    query_for_soma = f"""
            SELECT DISTINCT sp.id as specimen_id, 'null', layert.name as path_type, poly.path, sc.resolution, 'null', 'null'
            FROM specimens sp
            JOIN biospecimen_polygons AS bsp ON bsp.biospecimen_id=sp.id
            JOIN avg_graphic_objects poly ON poly.id=bsp.polygon_id
            JOIN avg_graphic_objects layer ON layer.id=poly.parent_id
            JOIN avg_group_labels layert ON layert.id=layer.group_label_id
            AND layert.prevent_missing_polygon_structure=false
            JOIN sub_images AS si ON si.id=layer.sub_image_id
            AND si.failed=false
            JOIN images AS im ON im.id=si.image_id
            JOIN slides AS s ON s.id=im.slide_id
            JOIN scans AS sc ON sc.slide_id=s.id
            AND sc.superseded=false
            JOIN treatments t ON t.id = im.treatment_id AND t.id = 300080909 --Why?
            WHERE sp.id IN ({ids_str})
            ORDER BY sp.id
            """
    # all results returned as 'invalid_data' with only soma coords and resolution
    _, cell_data = get_data(query_for_soma)
    soma_centers = {k: cell_data[k]["soma_center"] for k in cell_specimen_ids
                   if cell_data[k]["soma_center"] is not None}
    resolution = cell_data[int(cell_specimen_ids[1])]["resolution"]
    return soma_centers, resolution

def layer_info_from_snap_polygons_output(output, resolution=1):
    layers = {}
    for polygon in output["polygons"]:
        layers[polygon['name']] = {'bounds': Polygon(resolution*np.array(polygon['path']))}
    for surface in output["surfaces"]:
        name = surface['name']
        path = list(resolution*np.array(surface['path']))
        if name=='pia':
            pia_path = path
        elif name=='wm':
            wm_path = path
        else:
            path = LineString(path)
            layer, side = name.split('_')
            layers[layer][f"{side}_surface"] = path
    return layers, pia_path, wm_path

def get_layer_depths(point, layer_polys, pia_path, wm_path, depth_interp, dx_interp, dy_interp, step_size=1.0, max_iter=1000):
    pia_path = ensure_linestring(pia_path)
    wm_path = ensure_linestring(wm_path)
    in_layer = [
        layer for layer in layer_polys if
        layer_polys[layer]['bounds'].intersects(Point(*point)) # checks for common boundary or interior
    ]

    if len(in_layer) == 0:
        raise LayerDepthError("Point not found in any layer")
    elif len(in_layer) == 1:
        start_layer = in_layer[0]
    else:
        raise LayerDepthError(f"Overlapping layers: {in_layer}")
    layer_poly = layer_polys[start_layer]

    _, pia_side_dist = ld.step_from_node(
        point, depth_interp, dx_interp, dy_interp, layer_poly['pia_surface'], step_size, max_iter
    )
    _, wm_side_dist = ld.step_from_node(
            point, depth_interp, dx_interp, dy_interp, layer_poly['wm_surface'], -step_size, max_iter
    )
    _, pia_distance = ld.step_from_node(
            point, depth_interp, dx_interp, dy_interp, pia_path, step_size, max_iter
    )
    _, wm_distance = ld.step_from_node(
            point, depth_interp, dx_interp, dy_interp, wm_path, -step_size, max_iter
    )
    # if layer in ['Layer2', 'Layer3']:
    # add normalized L2-3 depth for human cells
    
    pia_direction = np.array([dx_interp(point), dy_interp(point)])
    pia_direction /= np.linalg.norm(pia_direction)

    layer_thickness = wm_side_dist + pia_side_dist
    cortex_thickness = pia_distance + wm_distance
    out = {
        'position':point,
        'layer_depth': pia_side_dist,
        'layer_thickness': layer_thickness,
        'normalized_layer_depth': pia_side_dist/layer_thickness,
        'normalized_depth': pia_distance/cortex_thickness, #vs depth_interp(point)?
        'absolute_depth': pia_distance,
        'cortex_thickness': cortex_thickness,
        'wm_distance': wm_distance,
        'layer':start_layer,
        'pia_direction':pia_direction,
        }
    return out

def get_depths_slice(focal_plane_image_series_id, soma_centers, resolution=1, step_size=1.0, max_iter=1000):
    # soma_centers, resolution = get_cell_soma_data(cell_id_list)
    soma_centers = {cell: resolution*np.array(position) for cell, position in soma_centers.items()}

    parser = Parser(args=[], input_data=dict(
        focal_plane_image_series_id=focal_plane_image_series_id))
    parser.args.pop('log_level')
    output = run_snap_polygons(**parser.args)

    layers, pia_path, wm_path = layer_info_from_snap_polygons_output(output, resolution)

    depth_field, gradient_field = get_depth_gradient_field(
            pia_path,
            wm_path,
    )

    # nearest will get closest non-nan value on grid
    # needed since some points on boundary are outside
    interp_params = dict(method="nearest")

    depth_interp = ld.setup_interpolator(
        depth_field, None, **interp_params)
    dx_interp = ld.setup_interpolator(
        gradient_field, "dx", **interp_params)
    dy_interp = ld.setup_interpolator(
        gradient_field, "dy", **interp_params)


    outputs = {}
    for name, point in soma_centers.items():
        try:
            outputs[name] = get_layer_depths(point, layers, pia_path, wm_path, depth_interp, dx_interp, dy_interp,
                                            step_size, max_iter)
        except LayerDepthError as exc:
            logger.warning((f"Failure getting depth info for cell {name}: {exc}"))
        except Exception:
            logger.exception(f"Failure getting depth info for cell {name}")

    if len(outputs)>0:
        # add the mean location of successful cells only
        mean_soma_loc = np.stack([soma_centers[cell] for cell in outputs]).mean(axis=0)
        outputs['mean'] = get_layer_depths(mean_soma_loc, layers, pia_path, wm_path, depth_interp, dx_interp, dy_interp,
                                                step_size, max_iter)
            
    return outputs

