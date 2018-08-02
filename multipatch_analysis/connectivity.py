"""
Prototype code for analyzing connectivity and synaptic properties between cell classes.


"""

from __future__ import print_function, division

from collections import OrderedDict
from multipatch_analysis.database import database as db
from multipatch_analysis.connection_strength import ConnectionStrength, get_amps, get_baseline_amps
from multipatch_analysis.morphology import Morphology
from multipatch_analysis import constants


def measure_connectivity(pairs, cell_groups, cell_classes):
    """Given a list of cell pairs and a dict that groups cells together by class,
    return a structure that describes connectivity between cell classes.
    """
    results = OrderedDict()
    for pre_class, pre_group in cell_groups.items():
        # inhibitory or excitatory class?
        pre_cre = cell_classes[pre_class].get('cre_type')
        pre_pyr = cell_classes[pre_class].get('pyramidal')
        class_is_exc = pre_cre == 'unknown' or pre_cre in constants.EXCITATORY_CRE_TYPES or pre_pyr is True

        for post_class, post_group in cell_groups.items():
            post_group = cell_groups[post_class]
            class_pairs = [p for p in pairs if p.pre_cell in pre_group and p.post_cell in post_group]
            probed_pairs = [p for p in class_pairs if pair_was_probed(p, class_is_exc)]
            connections_found = [p for p in probed_pairs if p.synapse]

            results[(pre_class, post_class)] = {
                'connections_found': connections_found,
                'pairs_probed': probed_pairs,
            }

            # if pre_class == 'L2/3 pyr' and post_class=='L2/3 pyr':
            #     raise Exception()
    
    return results


def classify_cells(cell_classes, cells=None, session=None):
    """Given cell class definitions and a list of cells, return a dict indicating which cells
    are members of each class.

    Parameters
    ----------
    cell_classes : dict
        Dict of {class_name: class_criteria}, where each *class_criteria* value describes selection criteria for a cell class.
    cells : list | None
        List of Cell instances to be classified.
    session: Session | None
        If *cells* is not provided, then a database session may be given instead from which
        cells will be selected.
    """
    if cells is None:
        cells = session.query(db.Cell, db.Cell.cre_type, db.Cell.target_layer, Morphology.pyramidal).join(Morphology)
    cell_groups = OrderedDict([(class_name, set()) for class_name in cell_classes])
    for cell in cells:
        for class_name, cell_class in cell_classes.items():
            if cell_in_class(cell, cell_class):
                cell_groups[class_name].add(cell)
    return cell_groups


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


def cell_in_class(cell, cell_class):
    morpho = cell.morphology
    for k, v in cell_class.items():
        if hasattr(cell, k):
            if getattr(cell, k) != v:
                return False
        elif hasattr(morpho, k):
            if getattr(morpho, k) != v:
                return False
        else:
            raise Exception('Cannot use "%s" for cell typing; attribute not found on cell or cell.morphology' % k)
    return True


def cell_class_name(cre_type=None, target_layer=None, pyramidal=None):
    name = []
    if target_layer is not None:
        name.append('L' + target_layer)
    if pyramidal is not None:
        name.append('pyr')
    if cre_type is not None:
        name.append(cre_type)
    return ' '.join(name)


