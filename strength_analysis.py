"""
Big question: what's the best way to measure synaptic strength / connectivity?

"""


from __future__ import print_function, division

from collections import OrderedDict
from constants import EXCITATORY_CRE_TYPES, INHIBITORY_CRE_TYPES

import argparse
import sys
import pyqtgraph as pg
import os
import pickle
import pyqtgraph.multiprocess as mp
import numpy as np
import scipy.stats

from connection_detection import MultiPatchExperimentAnalyzer, EvokedResponseGroup
from synaptic_dynamics import DynamicsAnalyzer
from experiment_list import ExperimentList
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import TraceList
from neuroanalysis.filter import bessel_filter
from neuroanalysis.event_detection import exp_deconvolve


def measure_peak(trace, sign, baseline=(0e-3, 9e-3), response=(11e-3, 17e-3)):
    baseline = trace.time_slice(*baseline).data.mean()
    if sign == '+':
        peak = trace.time_slice(*response).data.max()
    else:
        peak = trace.time_slice(*response).data.min()
    return peak - baseline


def measure_sum(trace, sign, baseline=(0e-3, 9e-3), response=(12e-3, 17e-3)):
    baseline = trace.time_slice(*baseline).data.sum()
    peak = trace.time_slice(*response).data.sum()
    return peak - baseline
        






if __name__ == '__main__':
    app = pg.mkQApp()
    #pg.dbg()
    
    expt_index = sys.argv[1]
    pre_id, post_id = map(int, sys.argv[2:4])
    
    # Load experiment index
    cache_file = 'expts_cache.pkl'
    expts = ExperimentList(cache=cache_file)

    expt = expts[expt_index]
    
    sign = '-' if expt.cells[pre_id].cre_type in INHIBITORY_CRE_TYPES else '+'
    print("sign:", sign)
    #analyzer = MultiPatchExperimentAnalyzer(expt.data)
    #pulses = analyzer.get_evoked_responses(pre_id, post_id, clamp_mode='ic', pulse_ids=[0])
    
    analyzer = DynamicsAnalyzer(expt, pre_id, post_id, align_to='spike')
    
    # collect all first pulse responses
    #responses = analyzer.amp_group
    
    # collect all events
    responses = analyzer.all_events

    
    n_responses = len(responses)
    
    # do exponential deconvolution on all responses
    deconv = TraceList()
    grid1 = PlotGrid()
    grid1.set_shape(2, 1)
    grid1[0, 0].setLabels(left=('PSP', 'V'))
    grid1[1, 0].setLabels(bottom=('time', 's'))
    
    results = OrderedDict()
    raw_group = EvokedResponseGroup()
    deconv_group = EvokedResponseGroup()
    
    if len(responses) == 0:
        raise Exception("No data found for this synapse")


    def input_filter(trace):
        bsub = trace - np.median(trace.time_slice(0, 10e-3).data)
        filt = bessel_filter(bsub, 1000.)
        return filt
    
    def deconv_filter(trace):
        dec = exp_deconvolve(trace, 15e-3)
        baseline = np.median(dec.time_slice(trace.t0, trace.t0+10e-3).data)
        deconv = bessel_filter(dec-baseline, 300.)
        return deconv
        

    order = np.argsort([t.start_time for t in responses.responses])
    with pg.ProgressDialog("Measuring amplitudes...", maximum=n_responses) as dlg:
        for i in range(n_responses):
            r = responses.responses[order[i]]
            b = responses.baselines[order[i]]            
            r.t0 = 0
            b.t0 = 0

            add_to_avg = True
            stim_name = r.parent.parent.meta['stim_name']
            if '200Hz' in stim_name or '100Hz' in stim_name:
                print("skipped in average:", r.parent.parent)
                add_to_avg = False
            
            # lowpass raw data
            filt = input_filter(r)
            base_filt = input_filter(b)
            grid1[0, 0].plot(filt.time_values, filt.data, pen=(255, 255, 255, 100))
            if add_to_avg:
                raw_group.add(filt, None)
            
            results.setdefault('raw_peak', []).append((
                measure_peak(filt, sign),
                measure_peak(base_filt, sign)
            ))
            
            #results.setdefault('raw_sum', []).append((
                #measure_sum(filt, sign),
                #measure_sum(base_filt, sign)
            #))
            
            # deconvolve
            deconv = deconv_filter(r)
            base_deconv = deconv_filter(b)
            grid1[1, 0].plot(deconv.time_values, deconv.data, pen=(255, 255, 255, 100))
            if add_to_avg:
                deconv_group.add(deconv, None)
            
            results.setdefault('deconv_peak', []).append((
                measure_peak(deconv, sign),
                measure_peak(base_deconv, sign)
            ))
            
            #results.setdefault('deconv_sum', []).append((
                #measure_sum(deconv, sign),
                #measure_sum(base_deconv, sign)
            #))
            
            dlg += 1
            if dlg.wasCanceled():
                break
        
        
    
    grid1.show()
    raw_mean = raw_group.mean()
    grid1[0, 0].plot(raw_mean.time_values, raw_mean.data, pen={'color': 'g', 'width': 2}, shadowPen={'color': 'k', 'width': 3}, antialias=True)
    
    deconv_mean = deconv_group.mean()
    grid1[1, 0].plot(deconv_mean.time_values, deconv_mean.data, pen={'color': 'g', 'width': 2}, shadowPen={'color': 'k', 'width': 3}, antialias=True)

    
    plts = PlotGrid()
    plts.set_shape(1, len(results))
    
    for i,k in enumerate(results):
        amps = np.array([x[0] for x in results[k]])
        bases = np.array([x[1] for x in results[k]])
        x = np.linspace(0.0, 1.0, len(amps))
        
        amps = amps / bases.std()
        bases = bases / bases.std()
        
        ks_s, ks_p = scipy.stats.ks_2samp(amps, bases)
        print("%s%sks_2samp: %g p=%g" % (k, ' '*(30-len(k)), ks_s, ks_p))
        
        plt = plts[0, i]
        plt.plot(x, amps, pen=None, symbol='o', symbolBrush=(255, 255, 0, 150), symbolPen=None)
        plt.plot(x, bases, pen=None, symbol='o', symbolBrush=(255, 0, 0, 150), symbolPen=None)
        plt.setTitle('%s<br>ks: %0.2g %0.2g' % (k, ks_s, ks_p))
        if i > 0:
            plt.hideAxis('left')
        
    plts.setXLink(plts[0,0])
    plts.setYLink(plts[0,0])
    plts[0, 0].setXRange(-0.2, 1.2)
    plts.show()

    

    
    