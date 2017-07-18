# *-* coding: utf-8 *-*

"""
Script comparing across multipatch experimental conditions.

"""

from __future__ import print_function, division

import argparse
import sys
import pyqtgraph as pg
import os
import pickle
import pyqtgraph.multiprocess as mp
import numpy as np

from synaptic_dynamics import DynamicsAnalyzer
from experiment_list import ExperimentList
from neuroanalysis.baseline import float_mode
from connection_detection import trace_mean
from neuroanalysis.data import Trace

parser = argparse.ArgumentParser()
parser.add_argument('--cre-type', type=str, help='Enter as pretype-posttype. If comparing connection types separate'
                    '' 'with ",". Ex pvalb-pavalb,pvalb-sst')
parser.add_argument('--calcium', action='store_true', default=False, dest='calcium',
                    help='cre-type must also be specified')
parser.add_argument('--age', type=str, help='Enter age ranges separated by ",". Ex 40-50,60-70.'
                    '' 'cre-type must also be specified')

args = parser.parse_args(sys.argv[1:])

cache_file = 'expts_cache.pkl'
all_expts = ExperimentList(cache=cache_file)
app = pg.mkQApp()
pg.dbg()

result_cache_file = 'connectivity_comparison_cache.pkl'
if os.path.exists(result_cache_file):
    try:
        result_cache = pickle.load(open(result_cache_file, 'rb'))
        print ('loaded cache')
    except:
        result_cache = {}
        sys.excepthook(*sys.exc_info())
        print ('Error loading cache')
else:
    result_cache = {}

def avg_first_pulse(expt, pre, post):
    key = (expt.nwb_file, pre, post)
    if key in result_cache:
        res = result_cache[key]
        if 'avg_est' not in res:
            return None, None
        avg_est = res['avg_est']
        avg_amp = Trace(data=res['data'], dt=res['dt'])
        return avg_est, avg_amp

    analyzer = DynamicsAnalyzer(expt, pre, post)
    avg_est, _, avg_amp, _ = analyzer.estimate_amplitude(plot=False)
    if avg_amp is None:
        result_cache[key] = {}
        ret = None, None
    else:
        result_cache[key] = {'avg_est': avg_est, 'data': avg_amp.data, 'dt': avg_amp.dt}
        ret = avg_est, avg_amp

    data = pickle.dumps(result_cache)
    open(result_cache_file, 'wb').write(data)
    print (key)
    return ret

def first_pulse_plot(expt_list, name=None):
    amp_plots = pg.plot(labels={'left': ('Vm', 'V')})
    amp_base_subtract = []
    for expt in expt_list:
        for pre, post in expt.connections:
            if expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]:
                avg_est, avg_amp = avg_first_pulse(expt, pre, post)
                if avg_amp is not None:
                    avg_amp.t0 = 0
                    base = float_mode(avg_amp.data[:int(10e-3 / avg_amp.dt)])
                    amp_base_subtract.append(avg_amp.copy(data=avg_amp.data - base))
                    amp_plots.plot(avg_amp.time_values, avg_amp.data - base)
                    app.processEvents()
    grand_mean = trace_mean(amp_base_subtract)
    amp_plots.addLegend()
    amp_plots.plot(grand_mean.time_values, grand_mean.data, pen={'color': 'g', 'width': 3}, name=name)

if args.cre_type is not None and len(args.cre_type.split(',')) == 1:
    cre_type = args.cre_type.split('-')
    if args.calcium is True:
        expts = all_expts.select(cre_type=cre_type, calcium='high')
        legend = ("%s->%s, calcium = 2.0mM " % (cre_type[0], cre_type[1]))
        dist_plots = expts.distance_plot(cre_type[0], cre_type[1], color=(0, 10), name=legend)
        first_pulse_plot(expts, name=legend)
        expts = all_expts.select(calcium='low')
        legend = ("%s->%s, calcium = 1.3mM " % (cre_type[0], cre_type[1]))
        expts.distance_plot(cre_type[0], cre_type[1], plots=dist_plots, color=(5, 10), name=("%s->%s, calcium = 1.3mM " %(cre_type[0], cre_type[1])))
        first_pulse_plot(expts, name=legend)
    elif args.age is not None:
        ages = args.age.split(',')
        if len(ages) < 2:
            print("Please specify more than one age range")
        expts = all_expts.select(age=ages[0])
        plots = expts.distance_plot(cre_type[0], cre_type[1], color=(0, 10), name=("%s->%s, age = P%s " %(cre_type[0], cre_type[1], ages[0])))
        expts = all_expts.select(age=ages[1])
        expts.distance_plot(cre_type[0], cre_type[1], plots=plots, color=(5, 10), name=("%s->%s, age = P%s " %(cre_type[0], cre_type[1], ages[1])))
elif args.cre_type is None and (args.calcium is not None or args.age is not None):
    print('Error: in order to compare across conditions a single cre-type connection must be specified')
else:
    cre_types = args.cre_type.split(',')
    plots = None
    for i, type in enumerate(cre_types):
        cre_type = type.split('-')
        expts = all_expts.select(cre_type=cre_type)
        plots = expts.distance_plot(cre_type[0], cre_type[1], plots=plots, color=(i, len(cre_types)*1.3))