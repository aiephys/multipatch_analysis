# -*- coding: utf-8 -*-

"""
Prototype code for analyzing connectivity and synaptic properties between cell classes.


"""
from __future__ import print_function, division

from collections import OrderedDict
import numpy as np
import pyqtgraph as pg
import pandas as pd
from statsmodels.stats.proportion import proportion_confint
import multipatch_analysis.database as db
# from first_pulse_deconvolved_amps import get_deconvolved_first_pulse_amps
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.baseline import float_mode


thermal_colormap = pg.ColorMap(
                    [0, 0.3333, 0.6666, 1],
                    [(255, 255, 255, 255), (255, 220, 0, 255), (185, 0, 0, 255), (0, 0, 0, 255)],
            )

class FormattableNumber(float):
    
    @property
    def si_format(self):
        return pg.siFormat(self)
       
    # all of the below formats assume that the number is entered with no scaling and in the form (scale)(unit)
    @property
    def mV(self):
        if np.isnan(self):
            formatted_value = ''
        else:
            value = self*1e3
            formatted_value = ("%0.2f mV" % value)
        return formatted_value

    @property
    def uV(self):
        if np.isnan(self):
            formatted_value = ''
        else:
            value = self*1e6
            formatted_value = (u"%0.f μV" % value)
        return formatted_value

    @property
    def pA(self):
        if np.isnan(self):
            formatted_value = ''
        else:
            value = self*1e12
            formatted_value = ("%0.2f pA" % value)
        return formatted_value

    @property
    def ms(self):
        if np.isnan(self):
            formatted_value = ''
        else:
            value = self*1e3
            formatted_value = ("%0.2f ms" % value)
        return formatted_value


class ConnectivityAnalyzer(object):
    class SignalHandler(pg.QtCore.QObject):
        """Because we can't subclass from both QObject and QGraphicsRectItem at the same time
        """
        sigOutputChanged = pg.QtCore.Signal(object) #self

    def __init__(self):
        self.name = 'connectivity'
        self.results = None
        self.group_results = None
        self._signalHandler = ConnectivityAnalyzer.SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged

        self.fields = [
            ('None', {}),
            ('probed', {}),
            ('connected', {}),
            ('gap_junction', {}),
            ('distance', {'mode': 'range', 'units': 'm', 'defaults': {
                'Max': 200e-6,
                'colormap': pg.ColorMap(
                    [0, 0.25, 0.5, 0.75, 1.0],
                    [(255,255,100,255), (255,100,0,255), (0,0,100,255), (140,0,0,255), (80,0,80,255)],
            )}}),
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
            ('gap_junction_probability', {'mode': 'range', 'defaults': {
                'colormap': pg.ColorMap(
                   [0, 0.01, 0.03, 0.1, 0.3, 1.0],
                    [(0,0,100, 255), (80,0,80, 255), (140,0,0, 255), (255,100,0, 255), (255,255,100, 255), (255,255,255, 255)],
            )}}),
            ]

        self.summary_stat = {
            'conn_no_data': self.metric_summary,
            'probed': self.metric_summary,
            'connected': self.metric_summary,
            'gap_junction': [self.metric_summary, self.metric_conf],
            'connection_probability': [self.metric_summary, self.metric_conf],
            'gap_junction_probability': [self.metric_summary, self.metric_conf],
            'matrix_completeness': [self.metric_summary, self.metric_conf],
            'distance': [self.metric_summary, self.metric_conf],
            
        }

    def metric_summary(self, x): 
        if x.name == 'conn_no_data':
            return all(x)
        if x.name == 'distance':
            return np.nanmean(x)
        if x.name in ['connected', 'probed', 'gap_junction']:
            return sum(filter(None, x))
        else:
            p = x.apply(pd.Series)
            p1 = p.apply(sum)
            connected = float(p1[0])
            probed = float(p1[1])

            if x.name == 'matrix_completeness':
                probed_progress = probed/80
                connected_progress = connected/6
                return np.clip(np.where(probed_progress > connected_progress, probed_progress, connected_progress), 0, 1)
            elif x.name.endswith('probability'):
                if probed == 0.:
                    return float('nan')
                else:
                    return connected/probed       

    def metric_conf(self, x):
        if x.name.endswith('probability'):
            p = x.apply(pd.Series)
            p1 = p.apply(sum)
            connected = float(p1[0])
            probed = float(p1[1])
            return connection_probability_ci(connected, probed)
        if x.name == 'distance':
            return [-np.nanstd(x), np.nanstd(x)]
        else:
            return float('nan')
    
    def invalidate_output(self):
        self.results = None
        self.group_results = None

    def measure(self, pair_groups):
        """Given a list of cell pairs and a dict that groups cells together by class,
        return a structure that describes connectivity of each cell pair.
        """    
        if self.results is not None:
            return self.results

        results = OrderedDict()

        for key, class_pairs in pair_groups.items():
            pre_class, post_class = key
            for pair in class_pairs:
                no_data = False
                probed = pair_was_probed(pair, pre_class.is_excitatory)
                if probed is False:
                    no_data = True

                connected = pair.synapse if probed is True else False
                gap = pair.electrical if probed is True else False
                distance = pair.distance if probed is True else float('nan')
                

                results[pair] = {
                'conn_no_data': no_data,
                'pre_class': pre_class,
                'post_class': post_class,
                'probed': probed,
                'connected': connected,
                'gap_junction': gap,
                'distance': distance,
                'connection_probability': [int(connected) if connected is not None else 0, int(probed) if probed is not None else 0],
                'gap_junction_probability': [int(gap) if gap is not None else 0, int(probed) if probed is not None else 0],
                'matrix_completeness': [int(connected) if connected is not None else 0, int(probed) if probed is not None else 0],
                
                }

        self.results = pd.DataFrame.from_dict(results, orient='index')

        return self.results

    def group_result(self):
        if self.group_results is not None:
            return self.group_results

        self.group_results = self.results.groupby(['pre_class', 'post_class']).agg(self.summary_stat)
        return self.group_results

    def output_fields(self):

        return self.fields

    def print_element_info(self, pre_class, post_class, element, field_name):
        connections = element[element['connected'] == True].index.tolist()
        print ("Connection type: %s -> %s" % (pre_class, post_class))
        print ("Connected Pairs:")
        for connection in connections:
            print ("\t %s" % (connection))
        gap_junctions = element[element['gap_junction'] == True].index.tolist()
        print ("Gap Junctions:")
        for gap in gap_junctions:
            print ("\t %s" % (gap))
        probed_pairs = element[element['probed'] == True].index.tolist()
        print ("Probed Pairs:")
        for probed in probed_pairs:
            print ("\t %s" % (probed))
        
    def plot_element_data(self, pre_class, post_class, element, field_name, color='g', trace_plt=None):
        summary = element.agg(self.summary_stat)  
        val = summary[field_name]['metric_summary']
        line = pg.InfiniteLine(val, pen={'color': color, 'width': 2}, movable=False)
        scatter = None
        baseline_window = int(db.default_sample_rate * 5e-3)
        traces = []
        point_data = []
        connections = element[element['connected'] == True].index.tolist()
        for pair in connections:
            cs = pair.connection_strength
            trace = cs.ic_average_response
            if trace is not None:
                x_offset = cs.ic_fit_xoffset
                trace = format_trace(trace, baseline_window, x_offset, align='psp')
                trace_plt.plot(trace.time_values, trace.data)
                traces.append(trace)
        grand_trace = TraceList(traces).mean()
        name = ('%s->%s, n=%d' % (pre_class, post_class, len(traces)))
        trace_plt.plot(grand_trace.time_values, grand_trace.data, pen={'color': color, 'width': 3}, name=name)
        trace_plt.setXRange(0, 20e-3)
        trace_plt.setLabels(left=('', 'V'), bottom=('Time from stimulus', 's'))
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
        self.name = 'strength'
        self.results = None
        self.group_results = None
        self.pair_items = {}
        self._signalHandler = ConnectivityAnalyzer.SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged

        self.summary_stat = {
        'ic_fit_amp_all': [self.metric_summary, self.metric_conf],
        'ic_fit_xoffset_all': [self.metric_summary, self.metric_conf],
        'ic_fit_rise_time_all': [self.metric_summary, self.metric_conf],
        'ic_fit_decay_tau_all': [self.metric_summary, self.metric_conf],
        'vc_fit_amp_all': [self.metric_summary, self.metric_conf],
        'vc_fit_xoffset_all': [self.metric_summary, self.metric_conf],
        'vc_fit_rise_time_all': [self.metric_summary, self.metric_conf],
        'vc_fit_decay_tau_all': [self.metric_summary, self.metric_conf],
        'ic_amp_first_pulse': [self.metric_summary, self.metric_conf],
        'ic_latency_first_pulse': [self.metric_summary, self.metric_conf],
        'ic_rise_time_first_pulse': [self.metric_summary, self.metric_conf],
        'ic_decay_tau_first_pulse': [self.metric_summary, self.metric_conf],
        'vc_amp_first_pulse': [self.metric_summary, self.metric_conf],
        'vc_latency_first_pulse': [self.metric_summary, self.metric_conf],
        'vc_rise_time_first_pulse': [self.metric_summary, self.metric_conf],
        'vc_decay_tau_first_pulse': [self.metric_summary, self.metric_conf],
        'strength_no_data': self.metric_summary,
        }

        self.fields = [
            ('None',{}),
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
        ]
        

    def invalidate_output(self):
        self.results = None
        self.group_results = None
        self.pair_items = {}

    def measure(self, pair_groups):
        """Given a list of cell pairs and a dict that groups cells together by class,
        return a structure that describes strength and kinetics of each cell pair.
        """  
        if self.results is not None:
            return self.results

        results = OrderedDict()
        
        for key, class_pairs in pair_groups.items():
            pre_class, post_class = key

            for pair in class_pairs:
                if pair.synapse is not True:
                    no_data = True
                    cs = None
                    fpf = None
                elif pair.synapse is True:
                    no_data = False
                    cs = pair.connection_strength
                    fpf = pair.avg_first_pulse_fit

                results[pair] = {
                'strength_no_data': no_data,
                'pre_class': pre_class,
                'post_class': post_class,
                'ic_fit_amp_all': cs.ic_fit_amp if cs is not None else float('nan'),
                'ic_fit_xoffset_all': cs.ic_fit_xoffset if cs is not None else float('nan'),
                'ic_fit_rise_time_all': cs.ic_fit_rise_time if cs is not None else float('nan'),
                'ic_fit_decay_tau_all': cs.ic_fit_decay_tau if cs is not None else float('nan'),
                'vc_fit_amp_all': cs.vc_fit_amp if cs is not None else float('nan'),
                'vc_fit_xoffset_all': cs.vc_fit_xoffset if cs is not None else float('nan'),
                'vc_fit_rise_time_all': cs.vc_fit_rise_time if cs is not None else float('nan'),
                'vc_fit_decay_tau_all': cs.vc_fit_decay_tau if cs is not None else float('nan'),
                'ic_amp_first_pulse': fpf.ic_amp if fpf is not None else float('nan'),
                'ic_latency_first_pulse': fpf.ic_latency if fpf is not None else float('nan'),
                'ic_rise_time_first_pulse': fpf.ic_rise_time if fpf is not None else float('nan'),
                'ic_decay_tau_first_pulse': fpf.ic_decay_tau if fpf is not None else float('nan'),
                'vc_amp_first_pulse': fpf.vc_amp if fpf is not None else float('nan'),
                'vc_latency_first_pulse': fpf.vc_latency if fpf is not None else float('nan'),
                'vc_rise_time_first_pulse': fpf.vc_rise_time if fpf is not None else float('nan'),
                'vc_decay_tau_first_pulse': fpf.vc_decay_tau if fpf is not None else float('nan'),
                }

        self.results = pd.DataFrame.from_dict(results, orient='index')

        return self.results

    def group_result(self):
        if self.group_results is not None:
            return self.group_results

        self.group_results = self.results.groupby(['pre_class', 'post_class']).agg(self.summary_stat)
        return self.group_results

    def output_fields(self):

        return self.fields

    def metric_summary(self, x):
        if x.name == 'strength_no_data':
            return all(x)
        else:
            return np.nanmean(x)

    def metric_conf(self, x):
        return [-np.nanstd(x), np.nanstd(x)]

    def print_element_info(self, pre_class, post_class, element, field_name=None):
        if field_name is not None:
            fn = field_name.split('_all')[0] if field_name.endswith('all') else field_name.split('_first_pulse')[0]
            units = [field[1].get('units', '') for field in self.fields if field[0] == field_name][0] 
            print ("Connection type: %s -> %s" % (pre_class, post_class))    
            print ("\t Grand Average %s = %s" % (field_name, pg.siFormat(element[field_name].mean(), suffix=units)))
            print ("Connected Pairs:")
            no_qc_data = []
            for pair, value in element[field_name].iteritems():
                if pair.synapse is not True:
                    continue
                if np.isnan(value):
                    no_qc_data.append(pair)
                else:
                    print ("\t %s" % (pair))
                    print ("\t\t Average %s: %s" % (field_name, pg.siFormat(value, suffix=units)))
            print("\t No QC Data:")
            for pair in no_qc_data:
                print ("\t\t %s" % (pair))

    def plot_element_data(self, pre_class, post_class, element, field_name, color='g', trace_plt=None):
        fn = field_name.split('_all')[0] if field_name.endswith('all') else field_name.split('_first_pulse')[0]
        val = element[field_name].mean()
        line = pg.InfiniteLine(val, pen={'color': color, 'width': 2}, movable=False)
        scatter = None
        baseline_window = int(db.default_sample_rate * 5e-3)
        values = []
        traces = []
        point_data = []
        for pair, value in element[field_name].iteritems():
            if pair.synapse is not True:
                continue
            if np.isnan(value):
                continue
            if field_name.endswith('all'):
                cs = pair.connection_strength
                trace = cs.ic_average_response if field_name.startswith('ic') else cs.vc_average_response
                x_offset = cs.ic_fit_xoffset if field_name.startswith('ic') else cs.vc_fit_xoffset
            elif field_name.endswith('first_pulse'):
                fpf = pair.avg_first_pulse_fit
                if fpf is None:
                    continue
                trace = fpf.ic_avg_psp_data if field_name.startswith('ic') else fpf.vc_avg_psp_data
                x_offset = fpf.ic_latency if field_name.startswith('ic') else fpf.vc_latency
            if trace is None:
                continue
            values.append(value)
            trace = format_trace(trace, baseline_window, x_offset, align='psp')
            trace_item = trace_plt.plot(trace.time_values, trace.data)
            point_data.append(pair)
            trace_item.pair = pair
            trace_item.curve.setClickable(True)
            trace_item.sigClicked.connect(self.trace_plot_clicked)
            traces.append(trace)
            self.pair_items[pair.id] = [trace_item]
        y_values = pg.pseudoScatter(np.asarray(values, dtype=float), spacing=1)
        scatter = pg.ScatterPlotItem(symbol='o', brush=(color + (150,)), pen='w', size=12)
        scatter.setData(values, y_values + 10., data=point_data)
        for point in scatter.points():
            pair_id = point.data().id
            self.pair_items[pair_id].append(point)
        scatter.sigClicked.connect(self.scatter_plot_clicked)
        grand_trace = TraceList(traces).mean()
        name = ('%s->%s, n=%d' % (pre_class, post_class, len(traces)))
        trace_plt.plot(grand_trace.time_values, grand_trace.data, pen={'color': color, 'width': 3}, name=name)
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
        self.name = 'dynamics'
        self.results = None
        self.group_results = None
        self._signalHandler = ConnectivityAnalyzer.SignalHandler()
        self.sigOutputChanged = self._signalHandler.sigOutputChanged

        self.summary_stat = {
            'dynamics_no_data': self.metric_summary,
            'pulse_ratio_8_1_50hz': [self.metric_summary, self.metric_conf],
            'pulse_ratio_2_1_50hz': [self.metric_summary, self.metric_conf],
            'pulse_ratio_5_1_50hz': [self.metric_summary, self.metric_conf],
            'pulse_ratio_9_1_125ms': [self.metric_summary, self.metric_conf],
            'pulse_ratio_9_1_250ms': [self.metric_summary, self.metric_conf],
            'pulse_ratio_9_1_500ms': [self.metric_summary, self.metric_conf],
            'pulse_ratio_9_1_1000ms': [self.metric_summary, self.metric_conf],
            'pulse_ratio_9_1_2000ms': [self.metric_summary, self.metric_conf],
            'pulse_ratio_9_1_4000ms': [self.metric_summary, self.metric_conf],
            'pulse_ratio_8_1_10hz': [self.metric_summary, self.metric_conf],
            'pulse_ratio_8_1_20hz': [self.metric_summary, self.metric_conf],
            'pulse_ratio_8_1_100hz': [self.metric_summary, self.metric_conf],
            'pulse_ratio_8_1_200hz': [self.metric_summary, self.metric_conf],
        }
        
        self.fields = [
            ('None', {}),
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
            ('pulse_ratio_8_1_10hz', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
            ('pulse_ratio_8_1_20hz', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
            ('pulse_ratio_8_1_100hz', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
            ('pulse_ratio_8_1_200hz', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
             ('pulse_ratio_9_1_125ms', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
             ('pulse_ratio_9_1_250ms', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
             ('pulse_ratio_9_1_500ms', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
             ('pulse_ratio_9_1_1000ms', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
             ('pulse_ratio_9_1_2000ms', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
             ('pulse_ratio_9_1_4000ms', {'mode': 'range', 'defaults': {
                'Min': 0, 
                'Max': 2, 
                'colormap': pg.ColorMap(
                [0, 0.5, 1.0],
                [(0, 0, 255, 255), (255, 255, 255, 255), (255, 0, 0, 255)],
            )}}),
            ]

    def invalidate_output(self):
        self.results = None
        self.group_results = None

    def metric_summary(self, x):
        if x.name == 'dynamics_no_data':
            return all(x)
        else:
            return np.nanmean(x)

    def metric_conf(self, x):
        return [-np.nanstd(x), np.nanstd(x)]

    def measure(self, pair_groups):
        """Given a list of cell pairs and a dict that groups cells together by class,
        return a structure that describes dynamics of each cell pair.
        """  
        if self.results is not None:
            return self.results

        results = OrderedDict()
        for key, class_pairs in pair_groups.items():
            pre_class, post_class = key
            
            for pair in class_pairs:
                if pair.synapse is False:
                    no_data = True
                    dynamics = None
                elif pair.synapse is True:
                    no_data = False
                    dynamics = pair.dynamics

                results[pair] = {
                'dynamics_no_data': no_data,
                'pre_class': pre_class,
                'post_class': post_class,
                'pulse_ratio_8_1_50hz': dynamics.pulse_ratio_8_1_50hz if dynamics is not None else float('nan'),
                'pulse_ratio_2_1_50hz': dynamics.pulse_ratio_2_1_50hz if dynamics is not None else float('nan'),
                'pulse_ratio_5_1_50hz': dynamics.pulse_ratio_5_1_50hz if dynamics is not None else float('nan'),
                'pulse_ratio_9_1_125ms': dynamics.pulse_ratio_9_1_125ms if dynamics is not None else float('nan'),
                'pulse_ratio_9_1_250ms': dynamics.pulse_ratio_9_1_250ms if dynamics is not None else float('nan'),
                'pulse_ratio_9_1_500ms': dynamics.pulse_ratio_9_1_500ms if dynamics is not None else float('nan'),
                'pulse_ratio_9_1_1000ms': dynamics.pulse_ratio_9_1_1000ms if dynamics is not None else float('nan'),
                'pulse_ratio_9_1_2000ms': dynamics.pulse_ratio_9_1_2000ms if dynamics is not None else float('nan'),
                'pulse_ratio_9_1_4000ms': dynamics.pulse_ratio_9_1_4000ms if dynamics is not None else float('nan'),
                'pulse_ratio_8_1_10hz': dynamics.pulse_ratio_8_1_10hz if dynamics is not None else float('nan'),
                'pulse_ratio_8_1_20hz': dynamics.pulse_ratio_8_1_20hz if dynamics is not None else float('nan'),
                'pulse_ratio_8_1_100hz': dynamics.pulse_ratio_8_1_100hz if dynamics is not None else float('nan'),
                'pulse_ratio_8_1_200hz': dynamics.pulse_ratio_8_1_200hz if dynamics is not None else float('nan'),
                }

        
        self.results = pd.DataFrame.from_dict(results, orient='index')
        
        return self.results

    def group_result(self):
        if self.group_results is not None:
            return self.group_results

        self.group_results = self.results.groupby(['pre_class', 'post_class']).agg(self.summary_stat)
        return self.group_results

    def output_fields(self):
       
        return self.fields

    def print_element_info(self, pre_class, post_class, element, field_name=None):
        if field_name is not None:
            print ("Connection type: %s -> %s" % (pre_class, post_class))    
            print ("\t Grand Average %s = %s" % (field_name, element[field_name].mean()))
            print ("Connected Pairs:")
            no_qc_data = []
            for pair, value in element[field_name].iteritems():
                if pair.synapse is not True:
                    continue
                if np.isnan(value):
                    no_qc_data.append(pair)
                else:
                    print ("\t %s" % (pair))
                    print ("\t\t Average %s: %0.2f" % (field_name, value))
            print("\t No QC Data:")
            for pair in no_qc_data:
                print ("\t\t %s" % (pair))

    def plot_element_data(self, pre_class, post_class, element, field_name, color='g', trace_plt=None):
        trace_plt = None
        val = element[field_name].mean()
        line = pg.InfiniteLine(val, pen={'color': color, 'width': 2}, movable=False)
        scatter = None
        baseline_window = int(db.default_sample_rate * 5e-3)
        values = []
        traces = []
        point_data = []
        for pair, value in element[field_name].iteritems():
            if np.isnan(value):
                continue
            traces = []
            if trace_plt is not None:
                trace = cs.ic_average_response if field_name.startswith('ic') else cs.vc_average_response
                x_offset = cs.ic_fit_latency if field_name.startswith('ic') else cs.vc_fit_latency
                trace = format_trace(trace, baseline_window, x_offset, align='psp')
                trace_plt.plot(trace.time_values, trace.data)
                traces.append(trace)
            values.append(value)
            y_values = pg.pseudoScatter(np.asarray(values, dtype=float), spacing=1)
            scatter = pg.ScatterPlotItem(symbol='o', brush=(color + (150,)), pen='w', size=12)
            scatter.setData(values, y_values + 10.)
            if trace_plt is not None:
                grand_trace = TraceList(traces).mean()
                trace_plt.plot(grand_trace.time_values, grand_trace.data, pen={'color': color, 'width': 3})
                units = 'V' if field_name.startswith('ic') else 'A'
                trace_plt.setXRange(0, 20e-3)
                trace_plt.setLabels(left=('', units), bottom=('Time from stimulus', 's'))
        return line, scatter

    def summary(self, results, metric):
        print('')


def format_trace(trace, baseline_win, x_offset, align='spike'):
    # align can be to the pre-synaptic spike (default) or the onset of the PSP ('psp')
    baseline = float_mode(trace[0:baseline_win])
    trace = Trace(data=(trace-baseline), sample_rate=db.default_sample_rate)
    if align == 'psp':
        trace.t0 = -x_offset
    return trace

def get_all_output_fields(analyzer_list):
    data_fields = []
    # confidence_fields = []
    for analyzer in analyzer_list:
        data_fields.extend(analyzer.output_fields())
        # confidence_fields.extend(analyzer.output_fields()['show_confidence'])

    return data_fields

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
    pairs = pairs.join(db.Experiment).outerjoin(db.Slice) ## don't want to drop all pairs if there aren't slice entries
    pairs = pairs.outerjoin(db.ConnectionStrength)
    
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