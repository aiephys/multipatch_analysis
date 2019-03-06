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
    class SignalHandler(pg.QtCore.QObject):
        """Because we can't subclass from both QObject and QGraphicsRectItem at the same time
        """
        sigOutputChanged = pg.QtCore.Signal(object) #self

    def __init__(self):
        self.results = None
        self._signalHandler = ConnectivityAnalyzer.SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged

    def invalidate_output(self):
        self.results = None

    def measure(self, pair_groups):
        """Given a list of cell pairs and a dict that groups cells together by class,
        return a structure that describes connectivity between cell classes.
        """    
        self.results = OrderedDict()
        for key, class_pairs in pair_groups.items():
            pre_class, post_class = key
            no_data = False
            probed_pairs = [p for p in class_pairs if pair_was_probed(p, pre_class.is_excitatory)]
            connections_found = [p for p in probed_pairs if p.synapse]
            gap_junctions = [p for p in probed_pairs if p.electrical]

            n_connected = len(connections_found)
            n_probed = len(probed_pairs)
            probed_progress = n_probed / 80.
            connected_progress = n_connected / 6.
            total_progress = np.clip(np.where(probed_progress > connected_progress, probed_progress, connected_progress), 0, 1)
            conf_interval = connection_probability_ci(n_connected, n_probed)
            conn_prob = float('nan') if n_probed == 0 else n_connected / n_probed
            if n_probed == 0:
                no_data = True

            self.results[(pre_class, post_class)] = {
                'no_data': no_data,
                'n_probed': n_probed,
                'n_connected': n_connected,
                'connection_probability': conn_prob,
                'confidence_interval': conf_interval,
                'connected_pairs': connections_found,
                'probed_pairs': probed_pairs,
                'matrix_completeness': total_progress,
                'distance_distribution': float('nan'),
                'gap_junctions': gap_junctions,
                'None': None
            }
        
        return self.results

    def output_fields(self):
        cmap = pg.ColorMap(
            [0, 0.01, 0.03, 0.1, 0.3, 1.0],
            [(0,0,100, 255), (80,0,80, 255), (140,0,0, 255), (255,100,0, 255), (255,255,100, 255), (255,255,255, 255)],
        )

        fields = {'color_by': [
            ('n_probed', {}),
            ('n_connected', {}),
            ('connection_probability', {'mode': 'range'}),
            ('matrix_completeness', {'mode': 'range'}),
            ('distance_distribution', {'mode': 'range'}),
            ('gap_junctions', {'mode': 'range'}),
            ],
            'show_confidence': [
            'None',
            'confidence_interval',
            ],
        }

        defaults = {'color_by': 'connection_probability', 
            'text': '{n_connected}/{n_probed}', 
            'show_confidence': 'confidence_interval', 
            'colormap': cmap,
            'min': 0,
            'max': 1}

        return fields, defaults

    def print_element_info(self, pre_class, post_class):
        connections = self.results[(pre_class, post_class)]['connected_pairs']
        print ("Connection type: %s -> %s" % (pre_class, post_class))
        print ("Connected Pairs:")
        for connection in connections:
            print ("\t %s" % (connection))
        probed_pairs = self.results[(pre_class, post_class)]['probed_pairs']
        print ("Probed Pairs:")
        for probed in probed_pairs:
            print ("\t %s" % (probed))

    def summary(self, results):
        total_connected = 0
        total_probed = 0
        for connectivity in results.values():
            total_connected += connectivity['n_connected']
            total_probed += connectivity['n_probed']

        print ("Total connected / probed\t %d / %d" % (total_connected, total_probed))


class StrengthAnalyzer(object):
    class SignalHandler(pg.QtCore.QObject):
        """Because we can't subclass from both QObject and QGraphicsRectItem at the same time
        """
        sigOutputChanged = pg.QtCore.Signal(object) #self

    def __init__(self):
        self.results = None
        self._signalHandler = ConnectivityAnalyzer.SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged

    def invalidate_output(self):
        self.results = None

    def measure(self, pair_groups):
        self.results = OrderedDict()
        for key, class_pairs in pair_groups.items():
            pre_class, post_class = key
            no_data = False
            connections = [p for p in class_pairs if p.synapse]
            if len(connections) == 0:
                no_data = True
            connection_strength = [c.connection_strength for c in connections] if len(connections) > 0 else float('nan')
            ic_amps = filter(None, [cs.ic_amp_mean for cs in connection_strength]) if len(connections) > 0 else float('nan')
            vc_amps = filter(None, [cs.vc_amp_mean for cs in connection_strength]) if len(connections) > 0 else float('nan')
            ic_latencies = filter(None, [cs.ic_latency_mean for cs in connection_strength]) if len(connections) > 0 else float('nan')
            
            self.results[(pre_class, post_class)] = {
                'no_data': no_data,
                'connected_pairs': connections,
                'n_connections': len(connections),
                'ic_mean_amp': np.mean(ic_amps),
                'ic_amp_stdev': [-np.std(ic_amps), np.std(ic_amps)],
                'vc_mean_amp': np.mean(vc_amps),
                'vc_amp_stdev': [-np.std(vc_amps), np.std(vc_amps)],
                'ic_mean_latency': np.mean(ic_latencies),
                'ic_latency_stdev': [-np.std(ic_latencies), np.std(ic_latencies)],
                'ic_mean_rise_time': float('nan'),
                'ic_rise_time_stdev': float('nan'),
                'ic_amp_cv': float('nan'),
                'None': None,
            }

        return self.results

    def output_fields(self):
        cmap = pg.ColorMap(
            [0, 0.5, 1.0],
            [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
        )

        fields = {'color_by': [
            ('n_connections', {}),
            ('ic_mean_amp', {'mode': 'range', 'units': 'V'}),
            ('vc_mean_amp', {'mode': 'range', 'units': 'A'}),
            ('ic_mean_latency', {'mode': 'range', 'units': 's'}),
            ('ic_mean_rise_time', {'mode': 'range', 'units': 's'}),
            ('ic_amp_cv', {'mode': 'range'}),
            ],
            'show_confidence': [
            'None',
            'ic_amp_stdev',
            'vc_amp_stdev',
            'ic_latency_stdev',
            'ic_rise_time_stdev'
            ],
        }

        defaults = {'color_by': 'ic_mean_amp', 
            'text': '{n_connections}', 
            'show_confidence': 'ic_amp_stdev', 
            'colormap': cmap,
            'min': -1e-3,
            'max': 1e-3}

        return fields, defaults

    def print_element_info(self, pre_class, post_class):
        connections = self.results[(pre_class, post_class)]['connected_pairs']
        print ("Connection type: %s -> %s" % (pre_class, post_class))
        print ("Connected Pairs:")
        for connection in connections:
            print ("\t %s" % (connection))
        
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