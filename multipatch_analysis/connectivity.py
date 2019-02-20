"""
Prototype code for analyzing connectivity and synaptic properties between cell classes.


"""
from __future__ import print_function, division

from collections import OrderedDict
import numpy as np
import pyqtgraph as pg
from statsmodels.stats.proportion import proportion_confint
from .database import database as db
from .connection_strength import ConnectionStrength
from .morphology import Morphology


class ConnectivityAnalyzer(object):
    def __init__(self, pair_groups):
        self.pair_groups = pair_groups

    def measure(self):
        """Given a list of cell pairs and a dict that groups cells together by class,
        return a structure that describes connectivity between cell classes.
        """    
        results = OrderedDict()
        for key, class_pairs in self.pair_groups.items():
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

    def display(self, pre_class, post_class, result):
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
        text = "%d/%d" % (result['n_connected'], result['n_probed'])
        output = set_display(connectivity, text, colormap, upper_ci=upper_ci, lower_ci=lower_ci)

        return output

    def summary(self, results):
        total_connected = 0
        total_probed = 0
        for connectivity in results.values():
            total_connected += connectivity['n_connected']
            total_probed += connectivity['n_probed']

        print ("Total connected / probed\t %d / %d" % (total_connected, total_probed))


class StrengthAnalyzer(object):
    def __init__(self, pair_groups):
        self.pair_groups = pair_groups

    def measure(self):
        results = OrderedDict()
        for key, class_pairs in self.pair_groups.items():
            pre_class, post_class = key

            connections = [p for p in class_pairs if p.synapse]
            if len(connections) > 0:
                connection_strength = [c.connection_strength for c in connections]
                results[(pre_class, post_class)] = {
                    'n_connections': len(connections),
                    'ic_mean_amp': [cs.ic_amp_mean for cs in connection_strength],
                    'ic_amp_stdev': [cs.ic_amp_stdev for cs in connection_strength],
                }
            else: 
                results[(pre_class, post_class)] = {'n_connections': 0}


        return results

    def display(self, pre_class, post_class, result):
        # this needs to be scaled by the range of amplitudes
        colormap = pg.ColorMap(
            [0, 0.5, 1.0],
            [(0, 0, 255), (255, 255, 255), (255, 0, 0)],
            )

        n_connections = result['n_connections']
        if n_connections > 0:
            amps = filter(None, result['ic_mean_amp'])
            grand_amp = np.mean(amps)
            text = ("%.2f mV" % (grand_amp*1e3))
            grand_stdev = np.std(amps)
            output = set_display(grand_amp, text, colormap, upper_ci=grand_stdev, lower_ci=(-grand_stdev))
        else:
            output = set_display(float('nan'), '', colormap)

        return output

    def summary(self, results):
        print('')

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
        if isinstance(acsf, str):
            pairs = pairs.filter(db.Experiment.acsf==acsf)
        else:
            pairs = pairs.filter(db.Experiment.acsf.in_(acsf))

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

def set_display(metric, text, colormap, upper_ci=None, lower_ci=None):
    if upper_ci is not None and lower_ci is not None:
        show_confidence = True
    else:
        show_confidence = False

    if show_confidence:
        output = {'bordercolor': 0.6}
        default_bgcolor = np.array([128., 128., 128.])
    else:
        output = {'bordercolor': 0.8}
        default_bgcolor = np.array([220., 220., 220.])
    
    if np.isnan(metric):
        output['bgcolor'] = tuple(default_bgcolor)
        output['fgcolor'] = 0.6
        output['text'] = ''
    else:
        # get color based on metric
        color = colormap.map(metric)
        
        # desaturate low confidence cells
        if show_confidence:
            confidence = (1.0 - (upper_ci - lower_ci)) ** 2
            color = color * confidence + default_bgcolor * (1.0 - confidence)
    
        # invert text color for dark background
        output['fgcolor'] = 'w' if sum(color[:3]) < 384 else 'k'
        output['text'] = text
        output['bgcolor'] = tuple(color)

    return output