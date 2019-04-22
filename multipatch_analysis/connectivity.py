"""
Prototype code for analyzing connectivity and synaptic properties between cell classes.


"""
from __future__ import print_function, division

from collections import OrderedDict
import numpy as np
import pyqtgraph as pg
from statsmodels.stats.proportion import proportion_confint
import multipatch_analysis.database as db
# from first_pulse_deconvolved_amps import get_deconvolved_first_pulse_amps
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.baseline import float_mode


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

    def plot_element_data(self, pre_class, post_class, field_name, data=None, color='g', trace_plt=None):
        val = self.results[(pre_class, post_class)][field_name]
        line = pg.InfiniteLine(val, pen={'color': color, 'width': 2}, movable=False)
        scatter = None
        baseline_window = int(db.default_sample_rate * 5e-3)
        if data is not None:
            values = []
            traces = []
            for pair in data:
                cs = pair.connection_strength
                values.append(getattr(cs, field_name))
                trace = cs.ic_average_response if field_name.startswith('ic') else cs.vc_average_response
                baseline = float_mode(trace[0:baseline_window])
                trace = Trace(data=(trace-baseline), sample_rate=db.default_sample_rate)
                trace_plt.plot(trace.time_values, trace.data)
                traces.append(trace)
            y_values = pg.pseudoScatter(np.asarray(values, dtype=float) + 1., spacing=0.3)
            scatter = pg.ScatterPlotItem(symbol='o', brush=(color + (150,)), pen='w', size=12)
            scatter.setData(values, y_values + 1.)
            grand_trace = TraceList(traces).mean()
            trace_plt.plot(grand_trace.time_values, grand_trace.data, pen={'color': color, 'width': 3})
            units = 'V' if field_name.startswith('ic') else 'A'
            trace_plt.setXRange(0, 20e-3)
            trace_plt.setLabels(left=('', units), bottom=('Time from stimulus', 's'))
        return line, scatter

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
        self.pair_items = {}
        self._signalHandler = ConnectivityAnalyzer.SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged

    def invalidate_output(self):
        self.results = None
        self.pair_items = {}

    def measure(self, pair_groups):
        self.results = OrderedDict()
        for key, class_pairs in pair_groups.items():
            pre_class, post_class = key
            no_data = False
            connections = [p for p in class_pairs if p.synapse]
            if len(connections) == 0:
                no_data = True
            # All pulses
            connection_strength = [c.connection_strength for c in connections] if len(connections) > 0 else float('nan')
            ic_amps_all = filter(None, [cs.ic_fit_amp for cs in connection_strength]) if len(connections) > 0 else float('nan')
            vc_amps_all = filter(None, [cs.vc_fit_amp for cs in connection_strength]) if len(connections) > 0 else float('nan')
            ic_xoffset_all = filter(None, [cs.ic_fit_xoffset for cs in connection_strength]) if len(connections) > 0 else float('nan')
            vc_xoffset_all = filter(None, [cs.vc_fit_xoffset for cs in connection_strength]) if len(connections) > 0 else float('nan')
            ic_rise_all = filter(None, [cs.ic_fit_rise_time for cs in connection_strength]) if len(connections) > 0 else float('nan')
            vc_rise_all = filter(None, [cs.vc_fit_rise_time for cs in connection_strength]) if len(connections) > 0 else float('nan')
            ic_decay_all = filter(None, [cs.ic_fit_decay_tau for cs in connection_strength]) if len(connections) > 0 else float('nan')
            vc_decay_all = filter(None, [cs.vc_fit_decay_tau for cs in connection_strength]) if len(connections) > 0 else float('nan')
            
            #First pulses
            first_pulse_fit = filter(None, [c.avg_first_pulse_fit for c in connections]) if len(connections) > 0 else float('nan')
            ic_amps_first = filter(None, [fpf.ic_amp for fpf in first_pulse_fit]) if len(connections) > 0 else float('nan')
            vc_amps_first = filter(None, [fpf.vc_amp for fpf in first_pulse_fit]) if len(connections) > 0 else float('nan')
            ic_latency_first = filter(None, [fpf.ic_latency for fpf in first_pulse_fit]) if len(connections) > 0 else float('nan')
            vc_latency_first = filter(None, [fpf.vc_latency for fpf in first_pulse_fit]) if len(connections) > 0 else float('nan')
            ic_rise_first = filter(None, [fpf.ic_rise_time for fpf in first_pulse_fit]) if len(connections) > 0 else float('nan')
            vc_rise_first = filter(None, [fpf.vc_rise_time for fpf in first_pulse_fit]) if len(connections) > 0 else float('nan')
            ic_decay_first = filter(None, [fpf.ic_decay_tau for fpf in first_pulse_fit]) if len(connections) > 0 else float('nan')
            vc_decay_first = filter(None, [fpf.vc_decay_tau for fpf in first_pulse_fit]) if len(connections) > 0 else float('nan')

            # cv = []
            # if len(connections) > 0:
            #     for connection in connections:
            #         fpda = get_deconvolved_first_pulse_amps(connection)
            #         if fpda is not None:
            #             cv.append(np.std(fpda)/np.mean(fpda))
            #         else:
            #             cv.append(float('nan'))
            # else:
            #     cv = float('nan')

            self.results[(pre_class, post_class)] = {
                'no_data': no_data,
                'connected_pairs': connections,
                'n_connections': len(connections),
                # 'cv_array': cv,
                ## all pulses
                'ic_fit_amp_all': np.mean(ic_amps_all),
                'ic_fit_xoffset_all': np.mean(ic_xoffset_all),
                'ic_fit_rise_time_all': np.mean(ic_rise_all),
                'ic_fit_decay_tau_all': np.mean(ic_decay_all),
                # 'ic_amp_cv': np.nanmedian(cv),
                'vc_fit_amp_all': np.mean(vc_amps_all),
                'vc_fit_xoffset_all': np.mean(vc_xoffset_all),
                'vc_fit_rise_time_all': np.mean(vc_rise_all),
                'vc_fit_decay_tau_all': np.mean(vc_decay_all),
                'ic_amp_stdev_all': [-np.std(ic_amps_all), np.std(ic_amps_all)],
                'ic_xoffset_stdev_all': [-np.std(ic_xoffset_all), np.std(ic_xoffset_all)],
                'ic_rise_time_stdev_all': [-np.std(ic_rise_all), np.std(ic_rise_all)],
                'ic_decay_tau_stdev_all': [-np.std(ic_decay_all), np.std(ic_decay_all)],
                'vc_amp_stdev_all': [-np.std(vc_amps_all), np.std(vc_amps_all)],
                'vc_xoffset_stdev_all': [-np.std(vc_xoffset_all), np.std(vc_xoffset_all)],
                'vc_rise_time_stdev_all': [-np.std(vc_rise_all), np.std(vc_rise_all)],
                'vc_decay_tau_stdev_all': [-np.std(vc_decay_all), np.std(vc_decay_all)],
                ## First pulse
                'ic_amp_first_pulse': np.nanmean(ic_amps_first),
                'ic_latency_first_pulse': np.nanmean(ic_latency_first),
                'ic_rise_time_first_pulse': np.nanmean(ic_rise_first),
                'ic_decau_tau_first_pulse': np.nanmean(ic_decay_first),
                'vc_amp_first_pulse': np.nanmean(vc_amps_first),
                'vc_latency_first_pulse': np.nanmean(vc_latency_first),
                'vc_rise_time_first_pulse': np.nanmean(vc_rise_first),
                'vc_decau_tau_first_pulse': np.nanmean(vc_decay_first),
                'ic_amp_stdev_first': [-np.std(ic_amps_first), np.std(ic_amps_first)],
                'ic_latency_stdev_first': [-np.std(ic_latency_first), np.std(ic_latency_first)],
                'ic_rise_time_stdev_first': [-np.std(ic_rise_first), np.std(ic_rise_first)],
                'ic_decay_rau_stdev_first': [-np.std(ic_decay_first), np.std(ic_decay_first)],
                'vc_amp_stdev_first': [-np.std(vc_amps_first), np.std(vc_amps_first)],
                'vc_latency_stdev_first': [-np.std(vc_latency_first), np.std(vc_latency_first)],
                'vc_rise_time_stdev_first': [-np.std(vc_rise_first), np.std(vc_rise_first)],
                'vc_decay_tau_stdev_first': [-np.std(vc_decay_first), np.std(vc_decay_first)],

                'None': None,
            }

        return self.results

    def output_fields(self):
       
        self.fields = {'color_by': [
            ('n_connections', {}),
            # all pulses
            ('ic_fit_amp_all', {'mode': 'range', 'units': 'V', 'defaults': {
                'Min': -1e-3, 
                'Max': 1e-3, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
            ('vc_fit_amp_all', {'mode': 'range', 'units': 'A', 'defaults': {
                'Min': -20e-12, 
                'Max': 20e-12, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(255, 0, 0, 255), (255, 255, 255, 255), (0, 0, 255, 255)],
            )}}),
            ('ic_fit_xoffset_all', {'mode': 'range', 'units': 's', 'defaults': {
                'Min': 0.5e-3, 
                'Max': 4e-3,
                'colormap': thermal_colormap,
            }}),
            ('vc_fit_xoffset_all', {'mode': 'range', 'units': 's', 'defaults': {
                'Min': 0.5e-3, 
                'Max': 4e-3,
                'colormap': thermal_colormap,
            }}),
            ('ic_fit_rise_time_all', {'mode': 'range', 'units': 's', 'defaults': {
                'Min': 1e-3, 
                'Max': 10e-3,
                'colormap': thermal_colormap,
            }}),
            ('vc_fit_rise_time_all', {'mode': 'range', 'units': 's', 'defaults': {
                'Min': 0.5e-3, 
                'Max': 5e-3,
                'colormap': thermal_colormap,
            }}),
            ('ic_fit_decay_tau_all', {'mode': 'range', 'units': 's', 'defaults': {
                'Max': 500e-3,
                'colormap': thermal_colormap,
            }}),
            ('vc_fit_decay_tau_all', {'mode': 'range', 'units': 's', 'defaults': {
                'Max': 20e-3,
                'colormap': thermal_colormap,
            }}),
            ## first pulse
            ('ic_amp_first_pulse', {'mode': 'range', 'units': 'V', 'defaults': {
                'Min': -1e-3, 
                'Max': 1e-3, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
            ('vc_amp_first_pulse', {'mode': 'range', 'units': 'A', 'defaults': {
                'Min': -20e-12, 
                'Max': 20e-12, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(255, 0, 0, 255), (255, 255, 255, 255), (0, 0, 255, 255)],
            )}}),
            ('ic_latency_first_pulse', {'mode': 'range', 'units': 's', 'defaults': {
                'Min': 0.5e-3, 
                'Max': 4e-3,
                'colormap': thermal_colormap,
            }}),
            ('vc_latency_first_pulse', {'mode': 'range', 'units': 's', 'defaults': {
                'Min': 0.5e-3, 
                'Max': 4e-3,
                'colormap': thermal_colormap,
            }}),
            ('ic_rise_time_first_pulse', {'mode': 'range', 'units': 's', 'defaults': {
                'Min': 1e-3, 
                'Max': 10e-3,
                'colormap': thermal_colormap,
            }}),
            ('vc_rise_time_first_pulse', {'mode': 'range', 'units': 's', 'defaults': {
                'Min': 0.5e-3, 
                'Max': 5e-3,
                'colormap': thermal_colormap,
            }}),
            ('ic_decay_tau_first_pulse', {'mode': 'range', 'units': 's', 'defaults': {
                'Max': 500e-3,
                'colormap': thermal_colormap,
            }}),
            ('vc_decay_tau_first_pulse', {'mode': 'range', 'units': 's', 'defaults': {
                'Max': 20e-3,
                'colormap': thermal_colormap,
            }}),

            # ('ic_amp_cv', {'mode': 'range', 'defaults': {
            #     'Min': -6,
            #     'Max': 6,
            #     'colormap': pg.ColorMap(
            #     [0, 0.5, 1.0],
            #     [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            # )}}),
            ],
            'show_confidence': [
            'None',
            'ic_amp_stdev_all',
            'vc_amp_stdev_all',
            'ic_xoffset_stdev_all',
            'vc_xoffset_stdev_all',
            'ic_rise_time_stdev_all',
            'vc_rise_time_stdev_all',
            'ic_decay_tau_stdev_all',
            'vc_decay_tau_stdev_all',
            'ic_amp_stdev_first',
            'vc_amp_stdev_first',
            'ic_latency_stdev_first',
            'vc_latency_stdev_first',
            'ic_rise_time_stdev_first',
            'vc_rise_time_stdev_first',
            'ic_decay_tau_stdev_first',
            'vc_decay_tau_stdev_first',
            ],
        }

        defaults = {'color_by': 'ic_fit_amp_all', 
            'text': '{n_connections}', 
            'show_confidence': 'ic_amp_stdev_all',
            'log': False,
        }

        return self.fields, defaults

    def print_element_info(self, pre_class, post_class, field_name=None):
        if field_name is not None:
            fn = field_name.split('_all')[0] if field_name.endswith('all') else field_name.split('_first_pulse')[0]
            units = [field[1].get('units', '') for field in self.fields['color_by'] if field[0] == field_name][0] 
            pairs = self.results[(pre_class, post_class)]['connected_pairs']
            result = self.results[(pre_class, post_class)][field_name]
            print ("Connection type: %s -> %s" % (pre_class, post_class))    
            print ("\t Grand Average %s = %s" % (field_name, pg.siFormat(result, suffix=units)))
            print ("Connected Pairs:")
            data = []
            for p, pair in enumerate(pairs):
                print ("\t %s" % (pair))
                # if metric == 'ic_amp_cv':
                #     value = self.results[(pre_class, post_class)]['cv_array'][c]
                #     pulses = 0
                # else:
                if field_name.endswith('all'):
                    cs = pair.connection_strength
                    value = getattr(cs, fn)
                    pulses = getattr(cs, fn.split('_')[0] + '_n_samples')
                elif field_name.endswith('first_pulse'):
                    fpf = pair.avg_first_pulse_fit
                    if fpf is not None:
                        value = getattr(fpf, fn)
                        pulses = len(fpf.ic_pulse_ids) if field_name.startswith('ic') else len(fpf.vc_pulse_ids)
                    else:
                        value = None
                        pulses = 0
                if value is not None:
                    print ("\t\t Average %s: %s, %d pulses" % (field_name, pg.siFormat(value, suffix=units), pulses))
                    data.append(pair)
                else:
                    print("\t\t No QC Data")

            return data

    def plot_element_data(self, pre_class, post_class, field_name, data=None, color='g', trace_plt=None):
        fn = field_name.split('_all')[0] if field_name.endswith('all') else field_name.split('_first_pulse')[0]
        val = self.results[(pre_class, post_class)][field_name]
        line = pg.InfiniteLine(val, pen={'color': color, 'width': 2}, movable=False)
        scatter = None
        baseline_window = int(db.default_sample_rate * 5e-3)
        if data is not None:
            values = []
            traces = []
            point_data = []
            for pair in data:
                if field_name.endswith('all'):
                    cs = pair.connection_strength
                    values.append(getattr(cs, fn))
                    trace = cs.ic_average_response if field_name.startswith('ic') else cs.vc_average_response
                elif field_name.endswith('first_pulse'):
                    fpf = pair.avg_first_pulse_fit
                    values.append(getattr(fpf, fn))
                    trace = fpf.ic_avg_psp_data if field_name.startswith('ic') else fpf.vc_avg_psp_data
                baseline = float_mode(trace[0:baseline_window])
                trace = Trace(data=(trace-baseline), sample_rate=db.default_sample_rate)
                trace_item = trace_plt.plot(trace.time_values, trace.data)
                point_data.append(pair)
                trace_item.pair = pair
                trace_item.curve.setClickable(True)
                trace_item.sigClicked.connect(self.trace_plot_clicked)
                traces.append(trace)
                self.pair_items[pair.id] = [trace_item]
            y_values = pg.pseudoScatter(np.asarray(values, dtype=float) + 1., spacing=0.3)
            scatter = pg.ScatterPlotItem(symbol='o', brush=(color + (150,)), pen='w', size=12)
            scatter.setData(values, y_values + 1., data=point_data)
            for point in scatter.points():
                pair_id = point.data().id
                self.pair_items[pair_id].append(point)
            scatter.sigClicked.connect(self.scatter_plot_clicked)
            grand_trace = TraceList(traces).mean()
            trace_plt.plot(grand_trace.time_values, grand_trace.data, pen={'color': color, 'width': 3})
            units = 'V' if field_name.startswith('ic') else 'A'
            trace_plt.setXRange(0, 20e-3)
            trace_plt.setLabels(left=('', units), bottom=('Time from stimulus', 's'))
        return line, scatter

    def scatter_plot_clicked(self, scatterplt, points):
        pair = points[0].data()
        self.select_pair(pair)

    def trace_plot_clicked(self, trace):
        pair = trace.pair
        self.select_pair(pair)

    def select_pair(self, pair):
        trace, point = self.pair_items[pair.id]
        point.setBrush(pg.mkBrush('y'))
        point.setSize(15)
        trace.setPen('y', width=2)
        print('Clicked:' '%s' % pair)

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
            pr_8_1_50 = []
            pr_2_1_50 = [] 
            pr_5_1_50 = []   
            dynamics = [c.dynamics for c in connections] if len(connections) > 0 else float('nan')
            for connection in connections:
                d = connection.dynamics
                pr_8_1_50.append(d.pulse_ratio_8_1_50hz if d is not None else float('nan'))
                pr_2_1_50.append(d.pulse_ratio_2_1_50hz if d is not None else float('nan'))
                pr_5_1_50.append(d.pulse_ratio_5_1_50hz if d is not None else float('nan'))


            self.results[(pre_class, post_class)] = {
            'no_data': no_data,
            'connected_pairs': connections,
            'n_connections': len(connections),
            'dynamics': dynamics,
            'pulse_ratio_8_1_50hz': np.nanmean(pr_8_1_50),
            'pulse_ratio_2_1_50hz': np.nanmean(pr_2_1_50),
            'pulse_ratio_5_1_50hz': np.nanmean(pr_5_1_50),
            'pulse_ratio_8_1_50hz_stdev': [-np.nanstd(pr_8_1_50), np.nanstd(pr_8_1_50)],
            'pulse_ratio_2_1_50hz_stdev': [-np.nanstd(pr_2_1_50), np.nanstd(pr_2_1_50)],
            'pulse_ratio_5_1_50hz_stdev': [-np.nanstd(pr_5_1_50), np.nanstd(pr_5_1_50)],
            'None': None,
            }
        
        return self.results

    def output_fields(self):
       
        self.fields = {'color_by': [
            ('pulse_ratio_8_1_50hz', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
            ('pulse_ratio_2_1_50hz', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
            ('pulse_ratio_5_1_50hz', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
            ],
            'show_confidence': [
            'pulse_ratio_8_1_50hz_stdev',
            'pulse_ratio_2_1_50hz_stdev',
            'pulse_ratio_5_1_50hz_stdev',
            'None',
            ]}

        defaults = {'color_by': 'pulse_ratio_8_1_50hz', 
            'text': '{n_connections}', 
            'show_confidence': 'None',
            'log': False,
        }

        return self.fields, defaults

    def print_element_info(self, pre_class, post_class, field_name=None):
        if field_name is not None:
            pairs = self.results[(pre_class, post_class)]['connected_pairs']
            result = self.results[(pre_class, post_class)][field_name]
            print ("Connection type: %s -> %s" % (pre_class, post_class))    
            print ("\t Grand Average %s = %s" % (field_name, result))
            print ("Connected Pairs:")
            data = []
            for pair in pairs:
                print ("\t %s" % (pair))
                d = pair.dynamics
                if d is not None:
                    value = getattr(d, field_name)
                    data.append(pair)
                    print ("\t\t Average %s: %s" % (field_name, value))
                else:
                    print("\t\t No QC Data")
        return data

    def plot_element_data(self, pre_class, post_class, field_name, data=None, color='g', trace_plt=None):
        val = self.results[(pre_class, post_class)][field_name]
        line = pg.InfiniteLine(val, pen={'color': color, 'width': 2}, movable=False)
        scatter = None
        baseline_window = int(db.default_sample_rate * 5e-3)
        if data is not None:
            values = []
            traces = []
            for pair in data:
                d = pair.dynamics
                values.append(getattr(d, field_name))
                if trace_plt is not None:
                    trace = cs.ic_average_response if field_name.startswith('ic') else cs.vc_average_response
                    baseline = float_mode(trace[0:baseline_window])
                    trace = Trace(data=(trace-baseline), sample_rate=db.default_sample_rate)
                    trace_plt.plot(trace.time_values, trace.data)
                    traces.append(trace)
            y_values = pg.pseudoScatter(np.asarray(values, dtype=float) + 1., spacing=0.3)
            scatter = pg.ScatterPlotItem(symbol='o', brush=(color + (150,)), pen='w', size=12)
            scatter.setData(values, y_values + 1.)
            if trace_plt is not None:
                grand_trace = TraceList(traces).mean()
                trace_plt.plot(grand_trace.time_values, grand_trace.data, pen={'color': color, 'width': 3})
                units = 'V' if field_name.startswith('ic') else 'A'
                trace_plt.setXRange(0, 20e-3)
                trace_plt.setLabels(left=('', units), bottom=('Time from stimulus', 's'))
        return line, scatter

    def summary(self, results, metric):
        print('')

def results_scatter(results, field_name, field, plt):
    vals = [result[field_name] for result in results.values() if np.isfinite(result[field_name])]
    y, x = np.histogram(vals, bins=np.linspace(min(vals), max(vals), 10))
    plt.plot(x, y, stepMode=True, fillLevel=0, brush=(255,255,255,150))
    units = field.get('units', '')
    plt.setLabels(left='Count', bottom=(field_name, units))

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
