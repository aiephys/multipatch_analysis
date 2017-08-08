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
import datetime
import re

from synaptic_dynamics import DynamicsAnalyzer
from experiment_list import ExperimentList
from neuroanalysis.baseline import float_mode
from connection_detection import trace_mean
from neuroanalysis.data import Trace
from scipy import stats
from neuroanalysis.ui.plot_grid import PlotGrid



def arg_to_date(arg):
    if arg is None:
        return None
    parts = re.split('\D+', arg)
    return datetime.date(*map(int, parts))

parser = argparse.ArgumentParser()
parser.add_argument('--cre-type', type=str, help='Enter as pretype-posttype. If comparing connection types separate'
                    '' 'with ",". Ex pvalb-pavalb,pvalb-sst')
parser.add_argument('--calcium', action='store_true', default=False, dest='calcium',
                    help='cre-type must also be specified')
parser.add_argument('--age', type=str, help='Enter age ranges separated by ",". Ex 40-50,60-70.'
                    '' 'cre-type must also be specified')
parser.add_argument('--start', type=arg_to_date)

args = parser.parse_args(sys.argv[1:])

cache_file = 'expts_cache.pkl'
all_expts = ExperimentList(cache=cache_file)
app = pg.mkQApp()

result_cache_file = 'synapse_comparison_cache.pkl'
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
            return None, None, None
        avg_est = res['avg_est']
        avg_amp = Trace(data=res['data'], dt=res['dt'])
        n_sweeps = res['n_sweeps']
        return avg_est, avg_amp, n_sweeps

    analyzer = DynamicsAnalyzer(expt, pre, post, align_to='pulse')
    avg_est, _, avg_amp, _, n_sweeps = analyzer.estimate_amplitude(plot=False)
    if n_sweeps == 0:
        result_cache[key] = {}
        ret = None, None, n_sweeps
    else:
        result_cache[key] = {'avg_est': avg_est, 'data': avg_amp.data, 'dt': avg_amp.dt, 'n_sweeps': n_sweeps}
        ret = avg_est, avg_amp, n_sweeps

    data = pickle.dumps(result_cache)
    open(result_cache_file, 'wb').write(data)
    print (key)
    return ret

def first_pulse_plot(expt_list, name=None):
    amp_plots = pg.plot()
    amp_plots.setLabels(left=('Vm', 'V'))
    amp_base_subtract = []
    avg_ests = []
    for expt in expt_list:
        for pre, post in expt.connections:
            if expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]:
                avg_est, avg_amp, n_sweeps = avg_first_pulse(expt, pre, post)
                if n_sweeps >= 10:
                    avg_amp.t0 = 0
                    avg_ests.append(avg_est)
                    base = float_mode(avg_amp.data[:int(10e-3 / avg_amp.dt)])
                    amp_base_subtract.append(avg_amp.copy(data=avg_amp.data - base))
                    amp_plots.plot(avg_amp.time_values, avg_amp.data - base)
                    app.processEvents()
    if len(amp_base_subtract) != 0:
        grand_mean = trace_mean(amp_base_subtract)
        grand_est = np.mean(np.array(avg_ests))
        amp_plots.addLegend()
        amp_plots.plot(grand_mean.time_values, grand_mean.data, pen={'color': 'g', 'width': 3}, name=name)
        amp_plots.addLine(y=grand_est, pen={'color': 'g'})
        return grand_mean, avg_ests, grand_est
    else:
        print ("No Traces")
        return None, avg_ests, None

def summary_plot(grand_mean, avg_est, grand_est, i, plot=None, color=None, name=None):
    if plot == None:
        grid = PlotGrid()
        grid.set_shape(1, 2)
        plot = (grid[0, 0], grid[0, 1])
        plot[0].grid = grid
        grid.show()
        plot[1].addLegend()
        plot[0].setLabels(left=('Vm', 'V'))
        plot[0].hideAxis('bottom')
        plot[1].setLabels(left=('Vm', 'V'))

    plot[1].plot(grand_mean.time_values, grand_mean.data, pen=color, name=name)
    dx = pg.pseudoScatter(np.array(avg_est).astype(float), 0.3, bidir=True)
    plot[0].plot((0.3 * dx / dx.max()) + i, avg_est, pen=None, symbol='x', symbolBrush=color,
                     symbolPen=None)
    plot[0].plot([i], [grand_est], pen=None, symbol='o', symbolBrush=color, symbolPen='w',
                     symbolSize=10)
    return plot

if args.cre_type is not None and len(args.cre_type.split(',')) == 1:
    cre_type = args.cre_type.split('-')
    if args.calcium is True:
        expts = all_expts.select(cre_type=cre_type, calcium='high', start=args.start)
        legend = ("%s->%s, calcium = 2.0mM " % (cre_type[0], cre_type[1]))
        dist_plots = expts.distance_plot(cre_type[0], cre_type[1], color=(0, 10), name=legend)
        grand_mean, avg_est_high, grand_est = first_pulse_plot(expts, name=legend)
        if grand_mean is not None:
            print(legend + 'Grand mean amplitude = %f' % grand_est)
            amp_plots = summary_plot(grand_mean, avg_est_high, grand_est, i=0, plot=None, color=(0, 10), name=legend)
        expts = all_expts.select(cre_type=cre_type, calcium='low', start=args.start)
        legend = ("%s->%s, calcium = 1.3mM " % (cre_type[0], cre_type[1]))
        expts.distance_plot(cre_type[0], cre_type[1], plots=dist_plots, color=(5, 10), name=legend)
        grand_mean, avg_est_low, grand_est = first_pulse_plot(expts, name=legend)
        if grand_mean is not None:
            print(legend + 'Grand mean amplitude = %f' % grand_est)
            amp_plots = summary_plot(grand_mean, avg_est_low, grand_est, i=1, plot=amp_plots, color=(5, 10), name=legend)
        ks = stats.ks_2samp(avg_est_high, avg_est_low)
        print('p = %f (KS test)' % ks.pvalue)
        #amp_plots[0].addLegend('p = %f (KS test)' % ks.pvalue)
    elif args.age is not None and len(args.age.split(',')) >= 2:
        ages = args.age.split(',')
        expts = all_expts.select(age=ages[0], start=args.start)
        legend = ("%s->%s, age = P%s " % (cre_type[0], cre_type[1], ages[0]))
        dist_plots = expts.distance_plot(cre_type[0], cre_type[1], color=(0, 10), name=legend)
        grand_mean, avg_est_age1, grand_est = first_pulse_plot(expts, name=legend)
        if grand_mean is not None:
            print(legend + 'Grand mean amplitude = %f' % grand_est)
            amp_plots = summary_plot(grand_mean, avg_est_age1, grand_est, i=0, plot=None, color=(0, 10), name=legend)
        expts = all_expts.select(age=ages[1], start=args.start)
        legend = ("%s->%s, age = P%s " % (cre_type[0], cre_type[1], ages[1]))
        expts.distance_plot(cre_type[0], cre_type[1], plots=dist_plots, color=(5, 10), name=legend)
        grand_mean, avg_est_age2, grand_est = first_pulse_plot(expts, name=legend)
        if grand_mean is not None:
            print(legend + 'Grand mean amplitude = %f' % grand_est)
            amp_plots = summary_plot(grand_mean, avg_est_age2, grand_est, i=1, plot=amp_plots, color=(5, 10), name=legend)
        ks = stats.ks_2samp(avg_est_age1, avg_est_age2)
        print('p = %f (KS test)' % ks.pvalue)
    elif args.cre_type is None and (args.calcium is not None or args.age is not None):
        print('Error: in order to compare across conditions a single cre-type connection must be specified')
    else:
        cre_types = args.cre_type.split(',')
        dist_plots = None
        amp_plots = None
        for i, type in enumerate(cre_types):
            cre_type = type.split('-')
            expts = all_expts.select(cre_type=cre_type, age=args.age, calcium='High', start=args.start)
            legend = ("%s->%s" % (cre_type[0], cre_type[1]))
            dist_plots = expts.distance_plot(cre_type[0], cre_type[1], plots=dist_plots, color=(i, len(cre_types)*1.3))
            grand_mean, avg_est, grand_est = first_pulse_plot(expts, name=legend)
            if grand_mean is not None:
                print(legend + 'Grand mean amplitude = %f' % grand_est)
                amp_plots = summary_plot(grand_mean, avg_est, grand_est, i=i, plot=amp_plots, color=(i, len(cre_types)*1.3), name=legend)