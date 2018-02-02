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

from multipatch_analysis.synaptic_dynamics import DynamicsAnalyzer
from multipatch_analysis.experiment_list import cached_experiments
from neuroanalysis.baseline import float_mode
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.filter import bessel_filter
from neuroanalysis.event_detection import exp_deconvolve
from scipy import stats
from neuroanalysis.ui.plot_grid import PlotGrid
from constants import INHIBITORY_CRE_TYPES, EXCITATORY_CRE_TYPES


def arg_to_date(arg):
    if arg is None:
        return None
    parts = re.split('\D+', arg)
    return datetime.date(*map(int, parts))

def load_cache(cache_file):
    if os.path.exists(cache_file):
        try:
            result_cache = pickle.load(open(cache_file, 'rb'))
            print ('loaded cache')
        except:
            result_cache = {}
            sys.excepthook(*sys.exc_info())
            print ('Error loading cache')
    else:
        result_cache = {}

    return result_cache

def responses(expt, pre, post):
    key = (expt.nwb_file, pre, post)
    result_cache = load_cache(result_cache_file)
    if key in result_cache:
        res = result_cache[key]
        if 'avg_est' not in res:
            return None, None, None
        avg_est = res['avg_est']
        avg_amp = Trace(data=res['data'], dt=res['dt'])
        n_sweeps = res['n_sweeps']
        return avg_est, avg_amp, n_sweeps

    analyzer = DynamicsAnalyzer(expt, pre, post, align_to='spike')
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

def first_pulse_plot(expt_list, name=None, summary_plot=None, color=None, scatter=0):
    amp_plots = pg.plot()
    amp_plots.setLabels(left=('Vm', 'V'))
    amp_base_subtract = []
    avg_ests = []
    for expt in expt_list:
        for pre, post in expt.connections:
            if expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]:
                avg_est, avg_amp, n_sweeps = responses(expt, pre, post)
                if expt.cells[pre].cre_type in EXCITATORY_CRE_TYPES and avg_est < 0:
                    continue
                elif expt.cells[pre].cre_type in INHIBITORY_CRE_TYPES and avg_est > 0:
                    continue
                if n_sweeps >= 10:
                    avg_amp.t0 = 0
                    avg_ests.append(avg_est)
                    base = float_mode(avg_amp.data[:int(10e-3 / avg_amp.dt)])
                    amp_base_subtract.append(avg_amp.copy(data=avg_amp.data - base))

                    current_connection_HS = post, pre
                    if len(expt.connections) > 1 and args.recip is True:
                        for i,x in enumerate(expt.connections):
                            if x == current_connection_HS:  # determine if a reciprocal connection
                                amp_plots.plot(avg_amp.time_values, avg_amp.data - base, pen={'color': 'r', 'width': 1})
                                break
                            elif x != current_connection_HS and i == len(expt.connections) - 1:  # reciprocal connection was not found
                                amp_plots.plot(avg_amp.time_values, avg_amp.data - base)
                    else:
                        amp_plots.plot(avg_amp.time_values, avg_amp.data - base)

                    app.processEvents()

    if len(amp_base_subtract) != 0:
        print(name + ' n = %d' % len(amp_base_subtract))
        grand_mean = TraceList(amp_base_subtract).mean()
        grand_est = np.mean(np.array(avg_ests))
        amp_plots.addLegend()
        amp_plots.plot(grand_mean.time_values, grand_mean.data, pen={'color': 'g', 'width': 3}, name=name)
        amp_plots.addLine(y=grand_est, pen={'color': 'g'})
        if grand_mean is not None:
            print(legend + ' Grand mean amplitude = %f' % grand_est)
            summary_plots = summary_plot_pulse(grand_mean, avg_ests, grand_est, labels=['Vm', 'V'], titles='Amplitude', i=scatter, plot=summary_plot, color=color, name=legend)
            return avg_ests, summary_plots
    else:
        print ("No Traces")
        return None, avg_ests, None, None

def train_response_plot(expt_list, name=None, summary_plots=[None, None], color=None):
    ind_base_subtract = []
    train_plots = pg.plot()
    train_plots.setLabels(left=('Vm', 'V'))
    tau =15e-3
    lp = 1000
    for expt in expt_list:
        for pre, post in expt.connections:
            if expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]:
                print ('Processing experiment: %s' % (expt.nwb_file))
                ind = []
                analyzer = DynamicsAnalyzer(expt, pre, post)
                train_responses = analyzer.train_responses
                for i, stim_params in enumerate(train_responses.keys()):
                     if stim_params[0] == 50:
                        if len(train_responses[stim_params][0]) != 0:
                            ind_group = train_responses[stim_params][0]
                            for j in range(len(ind_group)):
                                ind.append(ind_group.responses[j])
                if len(ind) > 5:
                    ind_avg = TraceList(ind).mean()
                    base = float_mode(ind_avg.data[:int(10e-3 / ind_avg.dt)])
                    ind_base_subtract.append(ind_avg.copy(data=ind_avg.data - base))
                    train_plots.plot(ind_avg.time_values, ind_avg.data - base)
                    app.processEvents()
    if len(ind_base_subtract) != 0:
        print (name + ' n = %d' % len(ind_base_subtract))
        ind_grand_mean = TraceList(ind_base_subtract).mean()
        ind_grand_mean_dec = bessel_filter(exp_deconvolve(ind_grand_mean, tau), lp)
        train_plots.addLegend()
        train_plots.plot(ind_grand_mean.time_values, ind_grand_mean.data, pen={'color': 'g', 'width': 3}, name=name)
        train_plots.plot(ind_grand_mean_dec.time_values, ind_grand_mean_dec.data, pen={'color': 'g', 'dash': [1,5,3,2]})
        if ind_grand_mean is not None:
            train_plots = summary_plot_train(ind_grand_mean, plot=summary_plots[0], color=color,
                                             name=(legend + ' 50 Hz induction'))
            train_plots2 = summary_plot_train(ind_grand_mean_dec, plot=summary_plots[1], color=color,
                                              name=(legend + ' 50 Hz induction'))
            return train_plots, train_plots2
    else:
        print ("No Traces")
        return None

def summary_plot_train(ind_grand_mean, plot=None, color=None, name=None):
    if plot is None:
        plot = pg.plot()
        plot.setLabels(left=('Vm', 'V'))
        plot.addLegend()

    plot.plot(ind_grand_mean.time_values, ind_grand_mean.data, pen=color, name=name)
    return plot

def summary_plot_pulse(feature_list, feature_mean, labels, titles, i, grand_trace=None, plot=None, color=None, name=None):
    if type(feature_list) is tuple:
        n_features = len(feature_list)
    else:
        n_features = 1
    if plot is None:
        plot = PlotGrid()
        plot.set_shape(n_features, 2)
        plot.show()
        for g in range(n_features):
            plot[g, 1].addLegend()
            plot[g, 1].setLabels(left=('Vm', 'V'))
            plot[g, 1].setLabels(bottom=('t', 's'))

    for feature in range(n_features):
        if n_features > 1:
            features = feature_list[feature]
            mean = feature_mean[feature]
            label = labels[feature]
            title = titles[feature]
        else:
            features = feature_list
            mean = feature_mean
            label = labels
            title = titles
        plot[feature, 0].setLabels(left=(label[0], label[1]))
        plot[feature, 0].hideAxis('bottom')
        plot[feature, 0].setTitle(title)
        if grand_trace is not None:
            plot[feature, 1].plot(grand_trace.time_values, grand_trace.data, pen=color, name=name)
        if len(features) > 1:
            dx = pg.pseudoScatter(np.array(features).astype(float), 0.3, bidir=True)
            bar = pg.BarGraphItem(x=[i], height=mean, width=0.7, brush='w', pen={'color': color, 'width': 2})
            plot[feature, 0].addItem(bar)
            sem = stats.sem(features)
            err = pg.ErrorBarItem(x=np.asarray([i]), y=np.asarray([mean]), height=sem, beam=0.3)
            plot[feature, 0].addItem(err)
            plot[feature, 0].plot((0.3 * dx / dx.max()) + i, features, pen=None, symbol='o', symbolSize=10, symbolPen='w',
                                symbolBrush=color)
        else:
            plot[feature, 0].plot([i], features, pen=None, symbol='o', symbolSize=10, symbolPen='w',
                                symbolBrush=color)

    return plot

def get_expts(all_expts, cre_type, calcium=None, age=None, start=None, dist_plot=None, color=None):
    expts = all_expts.select(cre_type=cre_type, calcium=calcium, age=age, start=start)
    if calcium is not None:
        if calcium == 'high':
            mM = '2.0 mM'
        elif calcium == 'low':
            mM = '1.3 mM'
    else:
        mM = ''
    if age is not None:
        age2 = 'P%s' % age
    else:
        age2 = 'All'
    legend = ("%s->%s, calcium = %s, age = %s" % (cre_type[0], cre_type[1], mM, age2))
    dist_plot = expts.distance_plot(cre_type[0], cre_type[1], plots=dist_plot, color=color, name=legend)
    return expts, legend, dist_plot

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cre-type', type=str, help='Enter as pretype-posttype. If comparing connection types separate'
                                                     '' 'with ",". Ex pvalb-pavalb,pvalb-sst')
    parser.add_argument('--calcium', action='store_true', default=False, dest='calcium',
                        help='cre-type must also be specified')
    parser.add_argument('--age', type=str, help='Enter age ranges separated by ",". Ex 40-50,60-70.'
                                                '' 'cre-type must also be specified')
    parser.add_argument('--recip', action='store_true', default=False, dest='recip',
                        help='Traces from reciprocal connections are red instead of gray')
    parser.add_argument('--start', type=arg_to_date)
    parser.add_argument('--trains', action='store_true', default=False, dest='trains',
                        help='optional to analyze 50Hz trains')

    args = parser.parse_args(sys.argv[1:])

    all_expts = cached_experiments()
    app = pg.mkQApp()

    result_cache_file = 'synapse_comparison_cache.pkl'
    if args.cre_type is not None:
        cre_types = args.cre_type.split(',')
        color = [(0, 10), (5, 10)]
        if args.calcium is True and len(cre_types) == 1:
            cre_type = args.cre_type.split('-')
            expts, legend, dist_plots = get_expts(all_expts, cre_type, calcium='high', age=args.age, start=args.start,
                                                  dist_plot=None, color=color[0])
            avg_est_high, summary_plots = first_pulse_plot(expts, name=legend, summary_plot=None, color=color[0], scatter=0)
            if args.trains is True:
                summary_train, summary_dec = train_response_plot(expts, name=(legend + ' 50 Hz induction'), summary_plots=[None,None],
                                                                 color=color[0])
            expts, legend, dist_plots = get_expts(all_expts, cre_type, calcium='low', age=args.age, start=args.start,
                                                    dist_plot=dist_plots, color=color[1])
            avg_est_low, summary_plots = first_pulse_plot(expts, name=legend, summary_plot=summary_plots, color=color[1], scatter=1)
            if args.trains is True:
                summary_train, summary_dec = train_response_plot(expts, name=(legend + ' 50 Hz induction'),
                                                                 summary_plots=[summary_train,summary_dec], color=color[1])
            ks = stats.ks_2samp(avg_est_high, avg_est_low)
            print('p = %f (KS test)' % ks.pvalue)
        elif args.age is not None and len(args.age.split(',')) >= 2 and len(cre_types) == 1:
            cre_type = args.cre_type.split('-')
            ages = args.age.split(',')
            expts, legend, dist_plots = get_expts(all_expts, cre_type, calcium=args.calcium, age=ages[0], start=args.start,
                                                  dist_plot=None, color=color[0])
            avg_est_age1, summary_plots = first_pulse_plot(expts, name=legend, summary_plot=None, color=color[0], scatter=0)
            if args.trains is True:
                summary_train, summary_dec = train_response_plot(expts, name=(legend + ' 50 Hz induction'), summary_plots=[None,None],
                                                                 color=color[0])
            expts, legend, dist_plots = get_expts(all_expts, cre_type, calcium=args.calcium, age=ages[1], start=args.start,
                                                  dist_plot=None, color=color[0])
            avg_est_age2, summary_plots = first_pulse_plot(expts, name=legend, summary_plot=summary_plots, color=color[1],
                                                           scatter=1)
            if args.trains is True:
                summary_train, summary_dec = train_response_plot(expts, name=(legend + ' 50 Hz induction'),
                                                                summary_plots=[summary_train, summary_dec], color=color[1])

            ks = stats.ks_2samp(avg_est_age1, avg_est_age2)
            print('p = %f (KS test)' % ks.pvalue)
        elif args.cre_type is None and (args.calcium is not None or args.age is not None):
            print('Error: in order to compare across conditions a single cre-type connection must be specified')
        else:
            dist_plots = None
            summary_plot = None
            train_plots = None
            train_plots2 = None
            for i, type in enumerate(cre_types):
                cre_type = type.split('-')
                expts, legend, dist_plots = get_expts(all_expts, cre_type, calcium='high', age=args.age, start=args.start,
                                                      dist_plot=dist_plots, color=(i, len(cre_types)*1.3))
                reciprocal_summary = expts.reciprocal(cre_type[0], cre_type[1])
                avg_est, summary_plot = first_pulse_plot(expts, name=legend, summary_plot=summary_plot, color=(i, len(cre_types)*1.3),
                                                         scatter=i)
                print(legend + ' %d/%d (%0.02f%%) uni-directional, %d/%d (%0.02f%%) reciprocal' % (
                  reciprocal_summary[tuple(cre_type)]['Uni-directional'], reciprocal_summary[tuple(cre_type)]['Total_connections'],
                  100 * reciprocal_summary[tuple(cre_type)]['Uni-directional']/reciprocal_summary[tuple(cre_type)]['Total_connections'],
                  reciprocal_summary[tuple(cre_type)]['Reciprocal'],
                  reciprocal_summary[tuple(cre_type)]['Total_connections'],
                  100 * reciprocal_summary[tuple(cre_type)]['Reciprocal']/reciprocal_summary[tuple(cre_type)]['Total_connections']))

                if args.trains is True:
                    train_plots, train_plots2 = train_response_plot(expts, name=(legend + ' 50 Hz induction'),
                                                                    summary_plots=[train_plots,train_plots2], color=(i, len(cre_types)*1.3))
