"""
Prototype code for analyzing connectivity and synaptic properties between cell classes.


"""

from __future__ import print_function, division

from collections import OrderedDict
from multipatch_analysis.database import database as db
from multipatch_analysis.connection_strength import ConnectionStrength, get_amps, get_baseline_amps
from multipatch_analysis import constants


def pair_was_probed(pair, excitatory):
    qc_field = 'n_%s_test_spikes' % ('ex' if excitatory else 'in')

    # arbitrary limit: we need at least N presynaptic spikes in order to consider
    # the pair "probed" for connection. Decreasing this value will decrease the number
    # of experiments included, but increase sensitivity for weaker connections
    return getattr(pair, qc_field) > 10


def cell_in_class(cell, cell_class):
    for k, v in cell_class.items():
        if getattr(cell, k, None) != v:
            return False
    return True


def cell_class_name(cre_type=None, target_layer=None):
    name = []
    if target_layer is not None:
        name.append('L' + target_layer)
    if cre_type is not None:
        name.append(cre_type)
    return ' '.join(name)


if __name__ == '__main__':

    session = db.Session()


    # 0. Define cell classes

    cell_classes = [
        {'cre_type': 'sim1'},
        {'cre_type': 'unknown', 'target_layer': '2/3'},
        {'cre_type': 'sst', 'target_layer': '2/3'},
        {'cre_type': 'pvalb', 'target_layer': '2/3'},
        {'cre_type': 'vip', 'target_layer': '2/3'},
        {'cre_type': 'sst', 'target_layer': '5'},
        {'cre_type': 'pvalb', 'target_layer': '5'},
        {'cre_type': 'vip', 'target_layer': '5'},
    ]
    cell_classes = OrderedDict([(cell_class_name(**cls), cls) for cls in cell_classes])


    # 1. Select pairs

    pre_cell = db.aliased(db.Cell, name='pre_cell')
    post_cell = db.aliased(db.Cell, name='post_cell')
    pairs = session.query(
        db.Pair, 
        pre_cell,
        post_cell,
        db.Experiment,
        db.Pair.synapse,
        ConnectionStrength.synapse_type,
    )
    pairs = pairs.join(pre_cell, pre_cell.id==db.Pair.pre_cell_id).join(post_cell, post_cell.id==db.Pair.post_cell_id)
    pairs = pairs.join(db.Experiment)
    pairs = pairs.join(ConnectionStrength)
    pairs = pairs.filter(db.Experiment.project_name=="mouse V1 coarse matrix")
    # calcium
    # age
    # egta
    
    pairs = pairs.all()


    # 2. Group cells by class

    cells = session.query(db.Cell, db.Cell.cre_type, db.Cell.target_layer)
    cell_groups = {}
    for cell in cells:
        for class_name, cell_class in cell_classes.items():
            if cell_in_class(cell, cell_class):
                cell_groups.setdefault(class_name, set()).add(cell.Cell.id)


    # 3. measure connectivity between groups

    for pre_class in cell_classes:
        # inhibitory or excitatory class?
        pre_cre = cell_classes[pre_class].get('cre_type')
        is_exc = pre_cre == 'unknown' or pre_cre in constants.EXCITATORY_CRE_TYPES
        pre_group = cell_groups[pre_class]

        for post_class in cell_classes:
            post_group = cell_groups[post_class]
            class_pairs = [p for p in pairs if p.pre_cell.id in pre_group and p.post_cell.id in post_group]
            probed_pairs = [p for p in class_pairs if pair_was_probed(p.Pair, is_exc)]
            connections_found = len([p for p in probed_pairs if p.synapse])
            connections_probed = len(probed_pairs)
            if connections_probed == 0:
                continue
            print("{pre_class:>20s} -> {post_class:20s} {connections_found} / {connections_probed}".format(
                pre_class=pre_class, 
                post_class=post_class, 
                connections_found=connections_found, 
                connections_probed=connections_probed,
            ))

