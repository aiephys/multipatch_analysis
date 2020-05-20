from __future__ import division
import os, pickle
from collections import OrderedDict
import numpy as np
import allensdk.core.swc as swc
import lims_utils
from layer_counts import layer_counts


def query_cells():
    q = """
    SELECT n.specimen_id, donors.full_genotype FROM neuron_reconstructions n 
    JOIN well_known_files f ON n.id = f.attachable_id
    JOIN specimens on specimens.id=n.specimen_id
    left join donors on specimens.donor_id=donors.id
    where 
        donors.full_genotype is not null
        and n.manual 
        and not n.superseded 
        and f.well_known_file_type_id = 303941301
    
    """
    cells = lims_utils.query(q, ())
    return cells


def soma_layer_from_morph(morph, layer_coords, soma_coords, layer_list):
    soma_x, soma_y = soma_coords["x"], soma_coords["y"]
    avg_x = soma_x.mean()
    avg_y = soma_y.mean()

    soma_node = morph.compartment_list_by_type(1)[0]
    delta_x = avg_x - soma_node["x"]
    delta_y = avg_y - soma_node["y"]

    soma_layer = layer_counts.point_in_layer(avg_x, avg_y, layer_coords, layer_list)
    return soma_layer


_morph_layer_density_cache = None
def morph_layer_density(specimen_id):
    global _morph_layer_density_cache
    cache_file = "morph_layer_density_cache.pkl"
    if _morph_layer_density_cache is None:
        if os.path.exists(cache_file):
            _morph_layer_density_cache = pickle.load(open(cache_file))
        else:
            _morph_layer_density_cache = {}
    
    key = specimen_id
    if key in _morph_layer_density_cache:
        return _morph_layer_density_cache[key]
        
    _morph_layer_density_cache[key] = _morph_layer_density(specimen_id)
    tmpfile = "cache.tmp"
    try:
        pickle.dump(_morph_layer_density_cache, open(tmpfile, 'w'))
        os.rename(tmpfile, cache_file)
    finally:
        if os.path.exists(tmpfile):
            os.remove(tmpfile)
            
    return _morph_layer_density_cache[key]


def _morph_layer_density(specimen_id, layer_list=None):
    layer_list = layer_list or ["1", "2/3", "4", "5", "6a", "6b"]

    df = layer_counts.query_lims_for_layers(specimen_id)
    try:
        soma_coords, pia_coords, wm_coords, layer_coords = layer_counts.get_coords(df, layer_list)
    except Exception:
        return None

    if None in [soma_coords, pia_coords, wm_coords, layer_coords]:
        return None

    ld = layer_counts.layer_depths(layer_coords, soma_coords, pia_coords, wm_coords)

    swc_filename, swc_path = lims_utils.get_swc_from_lims(specimen_id)

    morph = swc.read_swc(swc_path)

    soma_layer = soma_layer_from_morph(morph, layer_coords, soma_coords, layer_list)

    ax, dend, apic = layer_counts.compartment_counts_by_layer(morph, layer_coords, soma_coords, layer_list)
    
    all_dend = {k: dend[k] + apic[k] for k in dend}
    
    return {
        'soma_layer': soma_layer,
        'ax': ax,
        'dend': dend,
        'apic': apic,
        'all_dend': all_dend,
    }


def avg_density(cell_morpho_densities):
    avg = {}
    n_cells = len(cell_morpho_densities)
    for cell_morpho_density in cell_morpho_densities:
        for compartment_type, layers in cell_morpho_density.items():
            if compartment_type not in ['ax', 'dend', 'apic', 'all_dend']:
                continue
            avg.setdefault(compartment_type, {})
            for layer, density in layers.items():
                avg[compartment_type][layer] = avg[compartment_type].get(layer, 0) + density / n_cells
    
    return avg


def group_cells(cells):
    cells = cells[:]
    
    group_criteria = [
        ('2/3', 'cux2'),
        ('2/3', 'pvalb'),
        ('2/3', 'sst'),
        ('2/3', 'vip'),
        
        ('4', 'nr5a1'),
        ('4', 'pvalb'),
        ('4', 'sst'),
        ('4', 'vip'),
        
        ('5', 'sim1'),
        ('5', 'tlx3'),
        ('5', 'pvalb'),
        ('5', 'sst'),
        ('5', 'vip'),
        
        ('6a', 'ntsr1'),
        ('6a', 'pvalb'),
        ('6a', 'sst'),
        ('6a', 'vip'),
    ]
    groups = OrderedDict()
    
    for key in group_criteria:
        layer, cre = key
        groups[key] = []
        for cid, gtype in cells[:]:
            if not gtype.lower().startswith(cre):
                continue
            soma_layer = morph_layer_density(cid)['soma_layer']
            if layer != soma_layer:
                continue
            
            groups[key].append(cid)
            cells.remove((cid, gtype))
    
    return groups
            

def peters_rule_connectivity(dens1, dens2):
    pre = dens1['ax']
    post = dens2['all_dend']
    return sum([pre[k] * post[k] for k in pre])


def matrix_window(title="The Matrix"):
    """Create a new widget + viewbox for displaying matrix
    
    Returns
    -------
    w : GraphicsLayoutWidget
        The top-level widget created
    v : ViewBox
        The ViewBox embedded inside *w*
    """
    w = pg.GraphicsLayoutWidget()
    w.setRenderHints(w.renderHints() | pg.QtGui.QPainter.Antialiasing)
    w.setWindowTitle(title)
    v = w.addViewBox()
    v.setBackgroundColor('w')
    v.setAspectLocked()
    v.invertY()
    w.show()
    return w, v


if __name__ == '__main__':
    import pyqtgraph as pg
    pg.dbg()
    
    cells = query_cells()
    cells = [c for c in cells if morph_layer_density(c[0]) is not None]
    
    groups = group_cells(cells)
    
    avg_dens = {}
    for group_key, group_cells in groups.items():
        densities = [morph_layer_density(specimen_id) for specimen_id in group_cells]
        avg_dens[group_key] = avg_density(densities)

    for pre_group_key, pre_density in avg_dens.items():
        for post_group_key, post_density in avg_dens.items():
            conn = peters_rule_connectivity(pre_density, post_density)
            print(pre_group_key, post_group_key, conn)

    
    from multipatch_analysis.ui.graphics import MatrixItem
    
    w, v = matrix_window()

    colormap = pg.ColorMap(
        [0, 0.01, 0.03, 0.1, 0.3, 1.0],
        [(0,0,100), (80,0,80), (140,0,0), (255,100,0), (255,255,100), (255,255,255)],
    )

    shape = (len(groups),) * 2
    text = np.empty(shape, dtype=object)
    fgcolor = np.empty(shape, dtype=object)
    bgcolor = np.empty(shape, dtype=object)
    bordercolor = np.empty(shape, dtype=object)

    # call display function on every matrix element
    for i,pre_group_key in enumerate(groups):
        for j,post_group_key in enumerate(groups):
            pre_density = avg_dens[pre_group_key]
            post_density = avg_dens[post_group_key]
            conn = peters_rule_connectivity(pre_density, post_density)
            text[i, j] = ""
            fgcolor[i, j] = [0, 0, 0]
            bgcolor[i, j] = colormap.map(conn / 10e6)
            bordercolor[i, j] = [0, 0, 0]
            
    win, view = matrix_window(title="Peter's Rule Matrix")
    
    rows = list(groups.keys())
    matrix = MatrixItem(text=text, fgcolor=fgcolor, bgcolor=bgcolor, border_color=bordercolor,
                        rows=rows, cols=rows, size=50, header_color='k')
    view.addItem(matrix)
        
