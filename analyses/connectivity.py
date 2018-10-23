"""
Prototype code for analyzing connectivity and synaptic properties between cell classes.


"""

from __future__ import print_function, division

from collections import OrderedDict
import numpy as np
from multipatch_analysis.database import database as db
from multipatch_analysis.connectivity import query_pairs, measure_connectivity
from multipatch_analysis.connection_strength import ConnectionStrength, get_amps, get_baseline_amps
from multipatch_analysis.morphology import Morphology
from multipatch_analysis import constants
from multipatch_analysis.cell_class import CellClass, classify_cells, classify_pairs
from multipatch_analysis.ui.graphics import MatrixItem


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


def display_connectivity(pre_class, post_class, result, show_confidence=True):
    # Print results
    print("{pre_class:>20s} -> {post_class:20s} {connections_found:>5s} / {connections_probed}".format(
        pre_class=pre_class.name, 
        post_class=post_class.name, 
        connections_found=str(len(result['connected_pairs'])),
        connections_probed=len(result['probed_pairs']),
    ))

    # Pretty matrix results
    colormap = pg.ColorMap(
        [0, 0.01, 0.03, 0.1, 0.3, 1.0],
        [(0,0,100), (80,0,80), (140,0,0), (255,100,0), (255,255,100), (255,255,255)],
    )

    connectivity, lower_ci, upper_ci = result['connection_probability']

    if show_confidence:
        output = {'bordercolor': 0.6}
        default_bgcolor = np.array([128., 128., 128.])
    else:
        output = {'bordercolor': 0.8}
        default_bgcolor = np.array([220., 220., 220.])
    
    if np.isnan(connectivity):
        output['bgcolor'] = tuple(default_bgcolor)
        output['fgcolor'] = 0.6
        output['text'] = ''
    else:
        # get color based on connectivity
        color = colormap.map(connectivity)
        
        # desaturate low confidence cells
        if show_confidence:
            confidence = (1.0 - (upper_ci - lower_ci)) ** 2
            color = color * confidence + default_bgcolor * (1.0 - confidence)
        
        # invert text color for dark background
        output['fgcolor'] = 'w' if sum(color[:3]) < 384 else 'k'
        output['text'] = "%d/%d" % (result['n_connected'], result['n_probed'])
        output['bgcolor'] = tuple(color)

    return output


class MatrixAnalyzer(object):
    def __init__(self, cell_classes, pairs, analysis_func, display_func, title, session):
        self.matrix_window = None
        self.session = session

        # Group all cells by selected classes
        cell_groups = classify_cells(cell_classes, pairs=pairs)

        # Group pairs into (pre_class, post_class) groups
        pair_groups = classify_pairs(pairs, cell_groups)

        # analyze matrix elements
        results = analysis_func(pair_groups)

        shape = (len(cell_groups),) * 2
        text = np.empty(shape, dtype=object)
        fgcolor = np.empty(shape, dtype=object)
        bgcolor = np.empty(shape, dtype=object)
        bordercolor = np.empty(shape, dtype=object)

        # call display function on every matrix element
        for i,row in enumerate(cell_groups):
            for j,col in enumerate(cell_groups):
                output = display_func(row, col, results[(row, col)])
                text[i, j] = output['text']
                fgcolor[i, j] = output['fgcolor']
                bgcolor[i, j] = output['bgcolor']
                bordercolor[i, j] = output['bordercolor']
                
        # Force cell class descriptions down to tuples of 2 items
        # Kludgy, but works for now.
        rows = []
        for cell_class in cell_classes:
            tup = cell_class.as_tuple
            row = tup[:1]
            if len(tup) > 1:
                row = row + (' '.join(tup[1:]),)
            rows.append(row)

        win, view = matrix_window(title=title)
        self.matrix_window = win
        
        self.matrix = MatrixItem(text=text, fgcolor=fgcolor, bgcolor=bgcolor, border_color=bordercolor,
                            rows=rows, cols=rows, size=50, header_color='k')
        view.addItem(self.matrix)




if __name__ == '__main__':

    import pyqtgraph as pg
    pg.dbg()

    session = db.Session()
    
    # Define cell classes

    mouse_cell_classes = [
        # {'cre_type': 'unknown', 'pyramidal': True, 'target_layer': '2/3'},
        # {'cre_type': 'unknown', 'target_layer': '2/3'},
        # {'pyramidal': True, 'target_layer': '2/3'},
        {'pyramidal': True, 'target_layer': '2/3'},
        {'cre_type': 'sst', 'target_layer': '2/3'},
        {'cre_type': 'pvalb', 'target_layer': '2/3'},
        {'cre_type': 'vip', 'target_layer': '2/3'},
        {'cre_type': 'rorb', 'target_layer': '4'},
        {'cre_type': 'nr5a1', 'target_layer': '4'},
        {'cre_type': 'sst', 'target_layer': '4'},
        {'cre_type': 'pvalb', 'target_layer': '4'},
        {'cre_type': 'vip', 'target_layer': '4'},
        {'cre_type': 'sim1', 'target_layer': '5'},
        {'cre_type': 'tlx3', 'target_layer': '5'},
        {'cre_type': 'sst', 'target_layer': '5'},
        {'cre_type': 'pvalb', 'target_layer': '5'},
        {'cre_type': 'vip', 'target_layer': '5'},
        {'cre_type': 'ntsr1', 'target_layer': '6'},
        {'cre_type': 'sst', 'target_layer': '6'},
        {'cre_type': 'pvalb', 'target_layer': '6'},
        {'cre_type': 'vip', 'target_layer': '6'},
    ]

    human_cell_classes = [
        {'pyramidal': True, 'target_layer': '2'},
        {'pyramidal': False, 'target_layer': '2'},
        {'pyramidal': True, 'target_layer': '3'},
        {'pyramidal': False, 'target_layer': '3'},
        {'pyramidal': True, 'target_layer': '4'},
        {'pyramidal': False, 'target_layer': '4'},
        {'pyramidal': True, 'target_layer': '5'},
        {'pyramidal': False, 'target_layer': '5'},
        {'pyramidal': True, 'target_layer': '6'},
        {'pyramidal': False, 'target_layer': '6'},
    ]

    analyzers = []
    for cell_classes, project_names in [(mouse_cell_classes, ['mouse V1 coarse matrix', 'mouse V1 pre-production']), (human_cell_classes, ['human coarse matrix'])]:
        cell_classes = [CellClass(**c) for c in cell_classes]

        # Select pairs (todo: age, acsf, internal, temp, etc.)
        pairs = query_pairs(project_name=project_names, session=session).all()

        print("\n-------------------- %s ------------------\n" % ', '.join(project_names))
        maz = MatrixAnalyzer(cell_classes, pairs, analysis_func=measure_connectivity, display_func=display_connectivity, title=str(project_names), session=session)
        analyzers.append(maz)
