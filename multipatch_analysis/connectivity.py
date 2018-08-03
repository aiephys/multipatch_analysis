"""
Prototype code for analyzing connectivity and synaptic properties between cell classes.


"""

from __future__ import print_function, division

from collections import OrderedDict
from .database import database as db
from .connection_strength import ConnectionStrength
from .morphology import Morphology


def measure_connectivity(pairs, cell_groups):
    """Given a list of cell pairs and a dict that groups cells together by class,
    return a structure that describes connectivity between cell classes.
    """
    results = OrderedDict()
    for pre_class, pre_group in cell_groups.items():
        # inhibitory or excitatory class?

        for post_class, post_group in cell_groups.items():
            post_group = cell_groups[post_class]
            class_pairs = [p for p in pairs if p.pre_cell in pre_group and p.post_cell in post_group]
            probed_pairs = [p for p in class_pairs if pair_was_probed(p, pre_class.is_excitatory)]
            connections_found = [p for p in probed_pairs if p.synapse]

            results[(pre_class, post_class)] = {
                'connections_found': connections_found,
                'pairs_probed': probed_pairs,
            }

            # if pre_class == 'L2/3 pyr' and post_class=='L2/3 pyr':
            #     raise Exception()
    
    return results


def query_pairs(project_name, session):
    """Generate a query for selecting pairs from the database.

    Parameters
    ----------
    project_name : str
        Value to match from experiment.project_name (e.g. "mouse V1 coarse matrix" or "human coarse matrix")
    """
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
    # pairs = pairs.filter(db.Experiment.project_name=="mouse V1 coarse matrix")
    if project_name is not None:
        if isinstance(project_name, str):
            pairs = pairs.filter(db.Experiment.project_name==project_name)
        else:
            pairs = pairs.filter(db.Experiment.project_name.in_(project_name))
    # calcium
    # age
    # egta

    return pairs


def pair_was_probed(pair, excitatory):
    qc_field = 'n_%s_test_spikes' % ('ex' if excitatory else 'in')

    # arbitrary limit: we need at least N presynaptic spikes in order to consider
    # the pair "probed" for connection. Decreasing this value will decrease the number
    # of experiments included, but increase sensitivity for weaker connections
    return getattr(pair, qc_field) > 10

