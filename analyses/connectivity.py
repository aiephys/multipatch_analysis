"""
Prototype code for analyzing connectivity and synaptic properties between cell classes.


"""

from __future__ import print_function, division

from collections import OrderedDict
from multipatch_analysis.database import database as db
from multipatch_analysis.connectivity import query_pairs, classify_cells, measure_connectivity, cell_class_name
from multipatch_analysis.connection_strength import ConnectionStrength, get_amps, get_baseline_amps
from multipatch_analysis.morphology import Morphology
from multipatch_analysis import constants


if __name__ == '__main__':

    import pyqtgraph as pg
    pg.dbg()

    session = db.Session()


    # 0. Define cell classes

    mouse_cell_classes = [
        {'pyramidal': True, 'target_layer': '2/3'},
        {'cre_type': 'sst', 'target_layer': '2/3'},
        {'cre_type': 'pvalb', 'target_layer': '2/3'},
        {'cre_type': 'vip', 'target_layer': '2/3'},
        {'cre_type': 'rorb'},
        {'cre_type': 'nr5a1'},
        {'cre_type': 'sst', 'target_layer': '4'},
        {'cre_type': 'pvalb', 'target_layer': '4'},
        {'cre_type': 'vip', 'target_layer': '4'},
        {'cre_type': 'sim1'},
        {'cre_type': 'tlx3'},
        {'cre_type': 'sst', 'target_layer': '5'},
        {'cre_type': 'pvalb', 'target_layer': '5'},
        {'cre_type': 'vip', 'target_layer': '5'},
        {'cre_type': 'ntsr1'},
        {'cre_type': 'sst', 'target_layer': '6'},
        {'cre_type': 'pvalb', 'target_layer': '6'},
        {'cre_type': 'vip', 'target_layer': '6'},
    ]
    human_cell_classes = [
        {'pyramidal': True, 'target_layer': '2'},
        {'pyramidal': True, 'target_layer': '3'},
        {'pyramidal': True, 'target_layer': '4'},
        {'pyramidal': True, 'target_layer': '5'},
        {'pyramidal': True, 'target_layer': '6'},
    ]

    for cell_classes, project_name in [(mouse_cell_classes, 'mouse V1 coarse matrix'), (human_cell_classes, 'human coarse matrix')]:
        cell_classes = OrderedDict([(cell_class_name(**cls), cls) for cls in cell_classes])

        # 1. Select pairs
        records = query_pairs(project_name=project_name, session=session).all()
        pairs = [r[0] for r in records]

        # 2. Group all cells by selected classes
        cell_groups = classify_cells(cell_classes, session=session)

        # 3. measure connectivity between groups
        results = measure_connectivity(pairs, cell_groups, cell_classes)

        print("\n-------------------- %s ------------------\n" % project_name)
        for key, result in results.items():
            pre_class, post_class = key
            print("{pre_class:>20s} -> {post_class:20s} {connections_found:>5s} / {connections_probed}".format(
                pre_class=pre_class, 
                post_class=post_class, 
                connections_found=str(result['connections_found']),
                connections_probed=result['connections_probed'],
            ))
