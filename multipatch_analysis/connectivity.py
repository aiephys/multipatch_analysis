"""
Prototype code for analyzing connectivity and synaptic properties between cell classes.


"""
from __future__ import print_function, division

from collections import OrderedDict
import numpy as np
import pyqtgraph as pg
from statsmodels.stats.proportion import proportion_confint
import multipatch_analysis.database as db

thermal_colormap = pg.ColorMap(
                    [0, 0.3333, 0.6666, 1],
                    [(255, 255, 255, 255), (255, 220, 0, 255), (185, 0, 0, 255), (0, 0, 0, 255)],
            )

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
            n_gap_junctions = len(gap_junctions)
            n_probed = len(probed_pairs)
            probed_progress = n_probed / 80.
            connected_progress = n_connected / 6.
            total_progress = np.clip(np.where(probed_progress > connected_progress, probed_progress, connected_progress), 0, 1)
            conf_interval = connection_probability_ci(n_connected, n_probed)
            conn_prob = float('nan') if n_probed == 0 else n_connected / n_probed
            gap_prob = float('nan') if n_probed == 0 else n_gap_junctions / n_probed
            if n_probed == 0:
                no_data = True

            self.results[(pre_class, post_class)] = {
                'no_data': no_data,
                'n_probed': n_probed,
                'n_connected': n_connected,
                'connection_probability': conn_prob,
                'cp_confidence_interval': conf_interval,
                'connected_pairs': connections_found,
                'probed_pairs': probed_pairs,
                'matrix_completeness': total_progress,
                'distance_distribution': float('nan'),
                'gap_junctions': gap_junctions,
                'n_gap_junctions': n_gap_junctions,
                'gap_junction_probability': gap_prob,
                'None': None
            }
        
        return self.results

    def output_fields(self):

        self.fields = {'color_by': [
            ('n_probed', {}),
            ('n_connected', {}),
            ('n_gap_junctions', {}),
            ('connection_probability', {'mode': 'range', 'defaults': {
                'Operation': 'Add', 
                'colormap': pg.ColorMap(
                [0, 0.01, 0.03, 0.1, 0.3, 1.0],
                [(0,0,100, 255), (80,0,80, 255), (140,0,0, 255), (255,100,0, 255), (255,255,100, 255), (255,255,255, 255)],
            )}}),
            ('matrix_completeness', {'mode': 'range', 'defaults': {
                'colormap': pg.ColorMap(
                    [0, 0.25, 0.5, 0.75, 1.0],
                    [(0,0,100,255), (80,0,80,255), (140,0,0,255), (255,100,0,255), (255,255,100,255)],
            )}}),
            ('distance_distribution', {'mode': 'range'}),
            ('gap_junction_probability', {'mode': 'range', 'defaults': {
                'colormap': pg.ColorMap(
                   [0, 0.01, 0.03, 0.1, 0.3, 1.0],
                    [(0,0,100, 255), (80,0,80, 255), (140,0,0, 255), (255,100,0, 255), (255,255,100, 255), (255,255,255, 255)],
            )}}),
            ],
            'show_confidence': [
            'None',
            'cp_confidence_interval',
            ],
        }

        defaults = {'color_by': 'connection_probability', 
            'text': '{n_connected}/{n_probed}', 
            'show_confidence': 'cp_confidence_interval', 
            'log': True,
        }

        return self.fields, defaults

    def print_element_info(self, pre_class, post_class, metric):
        connections = self.results[(pre_class, post_class)]['connected_pairs']
        print ("Connection type: %s -> %s" % (pre_class, post_class))
        print ("Connected Pairs:")
        for connection in connections:
            print ("\t %s" % (connection))
        gap_junctions = self.results[(pre_class, post_class)]['gap_junctions']
        print ("Gap Junctions:")
        for gap in gap_junctions:
            print ("\t %s" % (gap))
        probed_pairs = self.results[(pre_class, post_class)]['probed_pairs']
        print ("Probed Pairs:")
        for probed in probed_pairs:
            print ("\t %s" % (probed))

    def summary(self, results, metric):
        if metric =='connection_probability':
            total_connected = 0
            total_probed = 0
            for connectivity in results.values():
                total_connected += connectivity['n_connected']
                total_probed += connectivity['n_probed']

            print ("Total connected / probed\t %d / %d" % (total_connected, total_probed))

        if metric == 'matrix_completeness':
            total_progress = 0
            for connectivity in results.values():
                total_progress += connectivity['matrix_completeness']
            n_elements = len([element for element in results.values() if element['no_data'] is False])

            print ("Total progress\t %0.1f%%, %d elements" % (100*total_progress/n_elements, n_elements))



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
            ic_amps = filter(None, [cs.ic_fit_amp for cs in connection_strength]) if len(connections) > 0 else float('nan')
            vc_amps = filter(None, [cs.vc_fit_amp for cs in connection_strength]) if len(connections) > 0 else float('nan')
            ic_xoffset = filter(None, [cs.ic_fit_xoffset for cs in connection_strength]) if len(connections) > 0 else float('nan')
            vc_xoffset = filter(None, [cs.vc_fit_xoffset for cs in connection_strength]) if len(connections) > 0 else float('nan')
            ic_rise = filter(None, [cs.ic_fit_rise_time for cs in connection_strength]) if len(connections) > 0 else float('nan')
            vc_rise = filter(None, [cs.vc_fit_rise_time for cs in connection_strength]) if len(connections) > 0 else float('nan')
            ic_decay = filter(None, [cs.ic_fit_decay_tau for cs in connection_strength]) if len(connections) > 0 else float('nan')
            vc_decay = filter(None, [cs.vc_fit_decay_tau for cs in connection_strength]) if len(connections) > 0 else float('nan')

            self.results[(pre_class, post_class)] = {
                'no_data': no_data,
                'connected_pairs': connections,
                'connection_strength': connection_strength,
                'n_connections': len(connections),
                'ic_fit_amp': np.mean(ic_amps),
                'ic_fit_xoffset': np.mean(ic_xoffset),
                'ic_fit_rise_time': np.mean(ic_rise),
                'ic_fit_decay_tau': np.mean(ic_decay),
                'vc_fit_amp': np.mean(vc_amps),
                'vc_fit_xoffset': np.mean(vc_xoffset),
                'vc_fit_rise_time': np.mean(vc_rise),
                'vc_fit_decay_tau': np.mean(vc_decay),
                'ic_amp_stdev': [-np.std(ic_amps), np.std(ic_amps)],
                'ic_xoffset_stdev': [-np.std(ic_xoffset), np.std(ic_xoffset)],
                'ic_rise_time_stdev': [-np.std(ic_rise), np.std(ic_rise)],
                'ic_decay_tau_stdev': [-np.std(ic_decay), np.std(ic_decay)],
                'vc_amp_stdev': [-np.std(vc_amps), np.std(vc_amps)],
                'vc_xoffset_stdev': [-np.std(vc_xoffset), np.std(vc_xoffset)],
                'vc_rise_time_stdev': [-np.std(vc_rise), np.std(vc_rise)],
                'vc_decay_tau_stdev': [-np.std(vc_decay), np.std(vc_decay)],
                'None': None,
            }

        return self.results

    def output_fields(self):
       
        self.fields = {'color_by': [
            ('n_connections', {}),
            ('ic_fit_amp', {'mode': 'range', 'units': 'V', 'defaults': {
                'Min': -1e-3, 
                'Max': 1e-3, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
            ('vc_fit_amp', {'mode': 'range', 'units': 'A', 'defaults': {
                'Min': -20e-12, 
                'Max': 20e-12, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(255, 0, 0, 255), (255, 255, 255, 255), (0, 0, 255, 255)],
            )}}),
            ('ic_fit_xoffset', {'mode': 'range', 'units': 's', 'defaults': {
                'Min': 0.5e-3, 
                'Max': 4e-3,
                'colormap': thermal_colormap,
            }}),
            ('vc_fit_xoffset', {'mode': 'range', 'units': 's', 'defaults': {
                'Min': 0.5e-3, 
                'Max': 4e-3,
                'colormap': thermal_colormap,
            }}),
            ('ic_fit_rise_time', {'mode': 'range', 'units': 's', 'defaults': {
                'Min': 1e-3, 
                'Max': 10e-3,
                'colormap': thermal_colormap,
            }}),
            ('vc_fit_rise_time', {'mode': 'range', 'units': 's', 'defaults': {
                'Min': 0.5e-3, 
                'Max': 5e-3,
                'colormap': thermal_colormap,
            }}),
            ('ic_fit_decay_tau', {'mode': 'range', 'units': 's', 'defaults': {
                'Max': 500e-3,
                'colormap': thermal_colormap,
            }}),
            ('vc_fit_decay_tau', {'mode': 'range', 'units': 's', 'defaults': {
                'Max': 20e-3,
                'colormap': thermal_colormap,
            }}),
            ],
            'show_confidence': [
            'None',
            'ic_amp_stdev',
            'vc_amp_stdev',
            'ic_latency_stdev',
            'ic_rise_time_stdev'
            ],
        }

        defaults = {'color_by': 'ic_fit_amp', 
            'text': '{n_connections}', 
            'show_confidence': 'ic_amp_stdev',
            'log': False,
        }

        return self.fields, defaults

    def print_element_info(self, pre_class, post_class, metric=None):
        if metric is not None:
            units = [field[1]['units'] for field in self.fields['color_by'] if field[0] == metric][0] 
            connections = self.results[(pre_class, post_class)]['connected_pairs']
            connection_strength = self.results[(pre_class, post_class)]['connection_strength']
            result = self.results[(pre_class, post_class)][metric]
            print ("Connection type: %s -> %s" % (pre_class, post_class))    
            print ("\t Grand Average %s = %s" % (metric, pg.siFormat(result, suffix=units)))
            print ("Connected Pairs:")
            values = []
            for connection in connections:
                print ("\t %s" % (connection))
                cs = [c for c in connection_strength if c.pair_id == connection.id][0]
                value = getattr(cs, metric)
                pulses = getattr(cs, metric.split('_')[0] + '_n_samples')
                if value is not None:
                    print ("\t\t Average %s: %s, %d pulses" % (metric, pg.siFormat(value, suffix=units), pulses))
                    values.append(value)
                else:
                    print("\t\t No QC Data")

            return values
        
    def summary(self, results, metric):
        print('')

class DynamicsAnalyzer(object):
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
            ind_50 = []
            ppr_50 = []    
            dynamics = [c.dynamics for c in connections] if len(connections) > 0 else float('nan')
            for connection in connections:
                d = connection.dynamics
                ind_50.append(d.pulse_ratio_8_1_50Hz if d is not None else float('nan'))
                ppr_50.append(d.pulse_ratio_2_1_50Hz if d is not None else float('nan'))

            self.results[(pre_class, post_class)] = {
            'no_data': no_data,
            'connected_pairs': connections,
            'n_connections': len(connections),
            'dynamics': dynamics,
            'pulse_ratio_8_1_50Hz': np.nanmean(ind_50),
            'pulse_ratio_2_1_50Hz': np.nanmean(ppr_50),
            'pulse_ratio_8_1_50Hz_stdev': [-np.nanstd(ind_50), np.nanstd(ind_50)],
            'pulse_ratio_2_1_50Hz_stdev': [-np.nanstd(ppr_50), np.nanstd(ppr_50)],
            'None': None,
            }
        
        return self.results

    def output_fields(self):
       
        self.fields = {'color_by': [
            ('pulse_ratio_8_1_50Hz', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
            ('pulse_ratio_2_1_50Hz', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
            ],
            'show_confidence': [
            'pulse_ratio_8_1_50Hz_stdev',
            'pulse_ratio_2_1_50Hz_stdev',
            'None',
            ]}

        defaults = {'color_by': 'pulse_ratio_8_1_50Hz', 
            'text': '{n_connections}', 
            'show_confidence': 'None',
            'log': False,
        }

        return self.fields, defaults

    def print_element_info(self, pre_class, post_class, metric=None):
        if metric is not None:
            connections = self.results[(pre_class, post_class)]['connected_pairs']
            result = self.results[(pre_class, post_class)][metric]
            print ("Connection type: %s -> %s" % (pre_class, post_class))    
            print ("\t Grand Average %s = %s" % (metric, result))
            print ("Connected Pairs:")
            values = []
            for connection in connections:
                print ("\t %s" % (connection))
                d = connection.dynamics
                if d is not None:
                    value = getattr(d, metric)
                    values.append(value)
                    print ("\t\t Average %s: %s" % (metric, value))
                else:
                    print("\t\t No QC Data")
        return values

    def summary(self, results, metric):
        print('')

def results_scatter(results, field_name, field, plt):
    vals = [result[field_name] for result in results.values() if np.isfinite(result[field_name])]
    y, x = np.histogram(vals, bins=np.linspace(min(vals), max(vals), 10))
    plt.plot(x, y, stepMode=True, fillLevel=0, brush=(255,255,255,150))
    units = field.get('units', '')
    plt.setLabels(left='Count', bottom=(field_name, units))

def add_element_to_scatter(results, pre_class, post_class, field_name, values=None, color='g'):
    val = results[(pre_class, post_class)][field_name]
    line = pg.InfiniteLine(val, pen={'color': color, 'width': 2}, movable=False)
    scatter = None
    if values is not None:
        y_values = pg.pseudoScatter(np.asarray(values, dtype=float) + 1., spacing=0.3)
        scatter = pg.ScatterPlotItem()
        scatter.setData(values, y_values + 1., symbol='o', brush=(color + (150,)), pen='w', size=12)
    return line, scatter

def connection_probability_ci(n_connected, n_probed):
    # make sure we are consistent about how we measure connectivity confidence intervals
    return proportion_confint(n_connected, n_probed, method='beta')

def query_pairs(project_name=None, acsf=None, age=None, species=None, distance=None, session=None, internal=None):
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
    pairs = pairs.join(db.ConnectionStrength)
    
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

    if internal is not None:
        if isinstance(internal, str):
            pairs = pairs.filter(db.Experiment.internal==internal)
        else:
            pairs = pairs.filter(db.Experiment.internal.in_(internal))

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