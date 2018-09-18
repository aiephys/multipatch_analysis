"""
Prototype code for analyzing connectivity and synaptic properties between cell classes.


"""
from __future__ import print_function, division

from collections import OrderedDict
import numpy as np
from statsmodels.stats.proportion import proportion_confint
from .database import database as db
from .connection_strength import ConnectionStrength
from .morphology import Morphology


def measure_connectivity(pair_groups):
    """Given a list of cell pairs and a dict that groups cells together by class,
    return a structure that describes connectivity between cell classes.
    """    
    results = OrderedDict()
    for key, class_pairs in pair_groups.items():
        pre_class, post_class = key
        
        probed_pairs = [p for p in class_pairs if pair_was_probed(p, pre_class.is_excitatory)]
        connections_found = [p for p in probed_pairs if p.synapse]

        n_connected = len(connections_found)
        n_probed = len(probed_pairs)
        conf_interval = connection_probability_ci(n_connected, n_probed)
        conn_prob = float('nan') if n_probed == 0 else n_connected / n_probed

        results[(pre_class, post_class)] = {
            'n_probed': n_probed,
            'n_connected': n_connected,
            'connection_probability': (conn_prob,) + conf_interval,
            'connected_pairs': connections_found,
            'probed_pairs': probed_pairs,
        }
    
    return results


def connection_probability_ci(n_connected, n_probed):
    # make sure we are consistent about how we measure connectivity confidence intervals
    return proportion_confint(n_connected, n_probed, method='beta')


def query_pairs(project_name=None, acsf=None, age=None, species=None, distance=None, session=None):
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
        # pre_cell,
        # post_cell,
        # db.Experiment,
        # db.Pair.synapse,
        # ConnectionStrength.synapse_type,
    )
    pairs = pairs.join(pre_cell, pre_cell.id==db.Pair.pre_cell_id).join(post_cell, post_cell.id==db.Pair.post_cell_id)
    pairs = pairs.join(db.Experiment).join(db.Slice)
    pairs = pairs.join(ConnectionStrength)
    
    if project_name is not None:
        if isinstance(project_name, str):
            pairs = pairs.filter(db.Experiment.project_name==project_name)
        else:
            pairs = pairs.filter(db.Experiment.project_name.in_(project_name))

    if acsf is not None:
        pairs = pairs.filter(db.Experiment.acsf==acsf)

    if age is not None:
        if age[0] is not None:
            pairs = pairs.filter(db.Slice.age>=age[0])
        if age[1] is not None:
            pairs = pairs.filter(db.Slice.age<=age[1])

    if distance is not None:
        if distance[0] is not None:
            pairs = pairs.filter(db.Pair.distance>=distance[0])
        if distance[1] is not None:
            pairs = pairs.filter(db.Pair.distance<=distance[1])

    if species is not None:
        pairs = pairs.filter(db.Slice.species==species)

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

