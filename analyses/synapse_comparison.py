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

from aisynphys.synaptic_dynamics import DynamicsAnalyzer
from aisynphys.experiment_list import cached_experiments
from neuroanalysis.baseline import float_mode
from neuroanalysis.data import TSeries, TSeriesList
from neuroanalysis.filter import bessel_filter
from neuroanalysis.event_detection import exp_deconvolve
from scipy import stats

from aisynphys.constants import INHIBITORY_CRE_TYPES, EXCITATORY_CRE_TYPES
from manuscript_figures import get_response, get_amplitude, response_filter, train_amp, write_cache
from aisynphys.fitting import fit_psp
from manuscript_figures import arg_to_date, load_cache, summary_plot_pulse, get_expts


def first_pulse_plot(expt_list, name=None, summary_plot=None, color=None, scatter=0, features=False):
    amp_plots = pg.plot()
    amp_plots.setLabels(left=('Vm', 'V'))
    grand_response = []
    avg_amps = {'amp': [], 'latency': [], 'rise': []}
    for expt in expt_list:
        if expt.connections is not None:
            for pre, post in expt.connections:
                if expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]:
                    all_responses, artifact = get_response(expt, pre, post, analysis_type='pulse')
                    if artifact > 0.03e-3:
                        continue
                    filtered_responses = response_filter(all_responses, freq_range=[0, 50], holding_range=[-68, -72], pulse=True)
                    n_sweeps = len(filtered_responses)
                    if n_sweeps >= 10:
                        avg_trace, avg_amp, amp_sign, _ = get_amplitude(filtered_responses)
                        if expt.cells[pre].cre_type in EXCITATORY_CRE_TYPES and avg_amp < 0:
                            continue
                        elif expt.cells[pre].cre_type in INHIBITORY_CRE_TYPES and avg_amp > 0:
                            continue
                        avg_trace.t0 = 0
                        avg_amps['amp'].append(avg_amp)
                        grand_response.append(avg_trace)
                        if features is True:
                            psp_fits = fit_psp(avg_trace, sign=amp_sign, yoffset=0, amp=avg_amp, method='leastsq',
                                               fit_kws={})
                            avg_amps['latency'].append(psp_fits.best_values['xoffset'] - 10e-3)
                            avg_amps['rise'].append(psp_fits.best_values['rise_time'])

                        current_connection_HS = post, pre
                        if len(expt.connections) > 1 and args.recip is True:
                            for i,x in enumerate(expt.connections):
                                if x == current_connection_HS:  # determine if a reciprocal connection
                                    amp_plots.plot(avg_trace.time_values, avg_trace.data, pen={'color': 'r', 'width': 1})
                                    break
                                elif x != current_connection_HS and i == len(expt.connections) - 1:  # reciprocal connection was not found
                                    amp_plots.plot(avg_trace.time_values, avg_trace.data)
                        else:
                            amp_plots.plot(avg_trace.time_values, avg_trace.data)

                        app.processEvents()

    if len(grand_response) != 0:
        print(name + ' n = %d' % len(grand_response))
        grand_mean = TSeriesList(grand_response).mean()
        grand_amp = np.mean(np.array(avg_amps['amp']))
        grand_amp_sem = stats.sem(np.array(avg_amps['amp']))
        amp_plots.addLegend()
        amp_plots.plot(grand_mean.time_values, grand_mean.data, pen={'color': 'g', 'width': 3}, name=name)
        amp_plots.addLine(y=grand_amp, pen={'color': 'g'})
        if grand_mean is not None:
            print(legend + ' Grand mean amplitude = %f +- %f' % (grand_amp, grand_amp_sem))
            if features is True:
                feature_list = (avg_amps['amp'], avg_amps['latency'], avg_amps['rise'])
                labels = (['Vm', 'V'], ['t', 's'], ['t', 's'])
                titles = ('Amplitude', 'Latency', 'Rise time')
            else:
                feature_list = [avg_amps['amp']]
                labels = (['Vm', 'V'])
                titles = 'Amplitude'
            summary_plots = summary_plot_pulse(feature_list[0], labels=labels, titles=titles, i=scatter,
                                               grand_trace=grand_mean, plot=summary_plot, color=color, name=legend)
            return avg_amps, summary_plots
    else:
        print ("No TSeries")
        return avg_amps, None

def train_response_plot(expt_list, name=None, summary_plots=[None, None], color=None):
    grand_train = [[], []]
    train_plots = pg.plot()
    train_plots.setLabels(left=('Vm', 'V'))
    tau =15e-3
    lp = 1000
    for expt in expt_list:
        for pre, post in expt.connections:
            if expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]:
                print ('Processing experiment: %s' % (expt.nwb_file))

                train_responses, artifact = get_response(expt, pre, post, analysis_type='train')
                if artifact > 0.03e-3:
                    continue

                train_filter = response_filter(train_responses['responses'], freq_range=[50, 50], train=0, delta_t=250)
                pulse_offsets = response_filter(train_responses['pulse_offsets'], freq_range=[50, 50], train=0, delta_t=250)

                if len(train_filter[0]) > 5:
                    ind_avg = TSeriesList(train_filter[0]).mean()
                    rec_avg = TSeriesList(train_filter[1]).mean()
                    rec_avg.t0 = 0.3
                    grand_train[0].append(ind_avg)
                    grand_train[1].append(rec_avg)
                    train_plots.plot(ind_avg.time_values, ind_avg.data)
                    train_plots.plot(rec_avg.time_values, rec_avg.data)
                    app.processEvents()
    if len(grand_train[0]) != 0:
        print (name + ' n = %d' % len(grand_train[0]))
        ind_grand_mean = TSeriesList(grand_train[0]).mean()
        rec_grand_mean = TSeriesList(grand_train[1]).mean()
        ind_grand_mean_dec = bessel_filter(exp_deconvolve(ind_grand_mean, tau), lp)
        train_plots.addLegend()
        train_plots.plot(ind_grand_mean.time_values, ind_grand_mean.data, pen={'color': 'g', 'width': 3}, name=name)
        train_plots.plot(rec_grand_mean.time_values, rec_grand_mean.data, pen={'color': 'g', 'width': 3}, name=name)
        train_amps = train_amp([grand_train[0], grand_train[1]], pulse_offsets, '+')
        if ind_grand_mean is not None:
            train_plots = summary_plot_train(ind_grand_mean, plot=summary_plots[0], color=color,
                                             name=(legend + ' 50 Hz induction'))
            train_plots = summary_plot_train(rec_grand_mean, plot=summary_plots[0], color=color)
            train_plots2 = summary_plot_train(ind_grand_mean_dec, plot=summary_plots[1], color=color,
                                              name=(legend + ' 50 Hz induction'))
            return train_plots, train_plots2, train_amps
    else:
        print ("No TSeries")
        return None

def summary_plot_train(ind_grand_mean, plot=None, color=None, name=None):
    if plot is None:
        plot = pg.plot()
        plot.setLabels(left=('Vm', 'V'))
        plot.addLegend()

    plot.plot(ind_grand_mean.time_values, ind_grand_mean.data, pen=color, name=name)
    return plot

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cre-type', type=str, help='Enter as pretype-posttype. If comparing connection types separate'
                                                     '' 'with ",". Ex pvalb-pavalb,pvalb-sst')
    parser.add_argument('--calcium', action='store_true', default=False, dest='calcium',
                        help='cre-type must also be specified')
    parser.add_argument('--age', type=str, help='Enter age ranges separated by ",". Ex 40-50,60-70.'
                                                '' 'cre-type must also be specified')
    parser.add_argument('--recip', action='store_true', default=False, dest='recip',
                        help='TSeries from reciprocal connections are red instead of gray')
    parser.add_argument('--start', type=arg_to_date)
    parser.add_argument('--stop', type=arg_to_date)
    parser.add_argument('--trains', action='store_true', default=False, dest='trains',
                        help='optional to analyze 50Hz trains')

    args = parser.parse_args(sys.argv[1:])

    all_expts = cached_experiments()
    app = pg.mkQApp()
    pg.dbg()
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')

    result_cache_file = 'synapse_comparison_cache.pkl'
    if args.cre_type is not None:
        cre_types = args.cre_type.split(',')
        color = [(255, 0, 0), (0, 0, 255)]
        if args.calcium is True and len(cre_types) == 1:
            cre_type = args.cre_type.split('-')
            expts, legend, dist_plots = get_expts(all_expts, cre_type, calcium='high', age=args.age, start=args.start,
                                                  stop=args.stop, dist_plot=None, color=color[0])
            avg_est_high, summary_plots = first_pulse_plot(expts, name=legend, summary_plot=None, color=color[0], scatter=0)
            if args.trains is True:
                summary_train, summary_dec = train_response_plot(expts, name=(legend + ' 50 Hz induction'), summary_plots=[None,None],
                                                                 color=color[0])
            expts, legend, dist_plots = get_expts(all_expts, cre_type, calcium='low', age=args.age, start=args.start,
                                                    stop=args.stop, dist_plot=dist_plots, color=color[1])
            avg_est_low, summary_plots = first_pulse_plot(expts, name=legend, summary_plot=summary_plots, color=color[1], scatter=1)
            if args.trains is True:
                summary_train, summary_dec = train_response_plot(expts, name=(legend + ' 50 Hz induction'),
                                                                 summary_plots=[summary_train,summary_dec], color=color[1])
            ks = stats.ks_2samp(avg_est_high['amp'], avg_est_low['amp'])
            print('p = %f (KS test)' % ks.pvalue)
        elif args.age is not None and len(args.age.split(',')) >= 2 and len(cre_types) == 1:
            cre_type = args.cre_type.split('-')
            ages = args.age.split(',')
            expts, legend, dist_plots = get_expts(all_expts, cre_type, calcium='high', age=ages[0], start=args.start,
                                                  stop=args.stop, dist_plot=None, color=color[0])
            avg_est_age1, summary_plots = first_pulse_plot(expts, name=legend, summary_plot=None, color=color[0], scatter=0,
                                                           features=True)
            if args.trains is True:
                summary_train, summary_dec = train_response_plot(expts, name=(legend + ' 50 Hz induction'), summary_plots=[None,None],
                                                                 color=color[0])
            expts, legend, dist_plots = get_expts(all_expts, cre_type, calcium='high', age=ages[1], start=args.start,
                                                  stop=args.stop, dist_plot=None, color=color[0])
            avg_est_age2, summary_plots = first_pulse_plot(expts, name=legend, summary_plot=summary_plots, color=color[1],
                                                           scatter=1, features=True)
            if args.trains is True:
                summary_train, summary_dec, summary_amps = train_response_plot(expts, name=(legend + ' 50 Hz induction'),
                                                                summary_plots=[summary_train, summary_dec], color=color[1])
                write_cache(summary_amps, 'age_train_amps.pkl')

            ks = stats.ks_2samp(avg_est_age1['amp'], avg_est_age2['amp'])
            print('p = %f (KS test, Amplitude)' % ks.pvalue)
            ks = stats.ks_2samp(avg_est_age1['latency'], avg_est_age2['latency'])
            print('p = %f (KS test, Latency)' % ks.pvalue)
            ks = stats.ks_2samp(avg_est_age1['rise'], avg_est_age2['rise'])
            print('p = %f (KS test, Rise time)' % ks.pvalue)
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
                                                      stop=args.stop, dist_plot=dist_plots, color=(i, len(cre_types)*1.3))
                reciprocal_summary = expts.reciprocal(cre_type[0], cre_type[1])
                avg_est, summary_plot = first_pulse_plot(expts, name=legend, summary_plot=summary_plot, color=(i, len(cre_types)*1.3),
                                                         scatter=i)
                if len(reciprocal_summary) > 0:
                    print(legend + ' %d/%d (%0.02f%%) uni-directional, %d/%d (%0.02f%%) reciprocal' % (
                      reciprocal_summary[tuple(cre_type)]['Uni-directional'], reciprocal_summary[tuple(cre_type)]['Total_connections'],
                      100 * reciprocal_summary[tuple(cre_type)]['Uni-directional']/reciprocal_summary[tuple(cre_type)]['Total_connections'],
                      reciprocal_summary[tuple(cre_type)]['Reciprocal'],
                      reciprocal_summary[tuple(cre_type)]['Total_connections'],
                      100 * reciprocal_summary[tuple(cre_type)]['Reciprocal']/reciprocal_summary[tuple(cre_type)]['Total_connections']))

                if args.trains is True:
                    train_plots, train_plots2 = train_response_plot(expts, name=(legend + ' 50 Hz induction'),
                                             summary_plots=[train_plots,train_plots2], color=(i, len(cre_types)*1.3))
