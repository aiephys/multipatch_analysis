"""
Prototype code for analyzing connectivity and synaptic properties between cell classes.


"""

from __future__ import print_function, division

from collections import OrderedDict
from multipatch_analysis.database import database as db
from multipatch_analysis.connectivity import query_pairs, measure_connectivity
from multipatch_analysis.connection_strength import ConnectionStrength, get_amps, get_baseline_amps
from multipatch_analysis.morphology import Morphology
from multipatch_analysis import constants
from multipatch_analysis.cell_class import CellClass, classify_cells


if __name__ == '__main__':

    import pyqtgraph as pg
    pg.dbg()

    session = db.Session()


    # 0. Define cell classes

    mouse_cell_classes = [
        {'cre_type': 'unknown', 'pyramidal': True, 'target_layer': '2/3'},
        {'cre_type': 'unknown', 'target_layer': '2/3'},
        {'pyramidal': True, 'target_layer': '2/3'},
        # {'cre_type': 'sst', 'target_layer': '2/3'},
        # {'cre_type': 'pvalb', 'target_layer': '2/3'},
        # {'cre_type': 'vip', 'target_layer': '2/3'},
        # {'cre_type': 'rorb'},
        # {'cre_type': 'nr5a1'},
        # {'cre_type': 'sst', 'target_layer': '4'},
        # {'cre_type': 'pvalb', 'target_layer': '4'},
        # {'cre_type': 'vip', 'target_layer': '4'},
        # {'cre_type': 'sim1'},
        # {'cre_type': 'tlx3'},
        # {'cre_type': 'sst', 'target_layer': '5'},
        # {'cre_type': 'pvalb', 'target_layer': '5'},
        # {'cre_type': 'vip', 'target_layer': '5'},
        # {'cre_type': 'ntsr1'},
        # {'cre_type': 'sst', 'target_layer': '6'},
        # {'cre_type': 'pvalb', 'target_layer': '6'},
        # {'cre_type': 'vip', 'target_layer': '6'},
    ]

    human_cell_classes = [
        {'pyramidal': True, 'target_layer': '2'},
        {'pyramidal': True, 'target_layer': '3'},
        {'pyramidal': True, 'target_layer': '4'},
        {'pyramidal': True, 'target_layer': '5'},
        {'pyramidal': True, 'target_layer': '6'},
    ]

    for cell_classes, project_names in [(mouse_cell_classes, ['mouse V1 coarse matrix', 'mouse V1 pre-production']), (human_cell_classes, ['human coarse matrix'])]:
        cell_classes = [CellClass(**c) for c in cell_classes]

        # 1. Select pairs (todo: age, acsf, internal, temp, etc.)
        records = query_pairs(project_name=project_names, session=session).all()
        pairs = [r[0] for r in records]
        cells = set([r[1] for r in records] + [r[2] for r in records])

        # 2. Group all cells by selected classes
        cell_groups = classify_cells(cell_classes, cells)

        # 3. measure connectivity between groups
        results = measure_connectivity(pairs, cell_groups)

        print("\n-------------------- %s ------------------\n" % ', '.join(project_names))
        for key, result in results.items():
            pre_class, post_class = key
            print("{pre_class:>20s} -> {post_class:20s} {connections_found:>5s} / {connections_probed}".format(
                pre_class=pre_class.name, 
                post_class=post_class.name, 
                connections_found=str(len(result['connections_found'])),
                connections_probed=len(result['pairs_probed']),
            ))
        break
