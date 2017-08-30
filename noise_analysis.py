"""
Big question: given that we only have so much time to record from a synapse,
what is the correct balance of repetitions versus number of different stimuli?
More repetitions gives us better confidence in our PSP amplitude estimates,
but more stimuli may be necessary to adequately sample the dynamic range of the
synapse.

Small question: What is the relationship between the number of repetitions and
the accuracy of PSP amplitude estimates?




"""


from __future__ import print_function, division

from collections import OrderedDict

import argparse
import sys
import pyqtgraph as pg
import os
import pickle
import pyqtgraph.multiprocess as mp
import numpy as np

from connection_detection import MultiPatchExperimentAnalyzer
from synaptic_dynamics import DynamicsAnalyzer
from experiment_list import ExperimentList
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import TraceList
from neuroanalysis.filter import bessel_filter
from neuroanalysis.event_detection import exp_deconvolve


if __name__ == '__main__':
    app = pg.mkQApp()
    #pg.dbg()
    
    expt_index, pre_id, post_id = map(int, sys.argv[1:4])
    
    # Load experiment index
    cache_file = 'expts_cache.pkl'
    expts = ExperimentList(cache=cache_file)

    expt = expts[expt_index]
    
    #analyzer = MultiPatchExperimentAnalyzer(expt.data)
    #pulses = analyzer.get_evoked_responses(pre_id, post_id, clamp_mode='ic', pulse_ids=[0])
    
    analyzer = DynamicsAnalyzer(expt, pre_id, post_id, align_to='spike')
    
    # collect all first pulse responses
    responses = analyzer.amp_group
    n_responses = len(responses)
    
    # do exponential deconvolution on all responses
    deconv = TraceList()
    grid1 = PlotGrid()
    grid1.set_shape(2, 1)
    for i in range(n_responses):
        r = responses.responses[i]
        grid1[0, 0].plot(r.time_values, r.data)
        
        dec = exp_deconvolve(r, 15e-3)
        baseline = np.median(dec.data[:100])
        r2 = bessel_filter(dec-baseline, 300.)
        grid1[1, 0].plot(r2.time_values, r2.data)
        
        deconv.append(r2)
    
    grid1.show()
    

    def measure_amp(trace, baseline=(6e-3, 8e-3), response=(12e-3, 20e-3)):
        baseline = trace.time_slice(*baseline).data.mean()
        peak = trace.time_slice(*response).data.max()
        return peak - baseline

    
    # Chop up responses into groups of varying size and plot the average
    # amplitude as measured from these chunks. We expect to see variance
    # decrease as a function of the number of responses being averaged together.
    x = []
    y_deconv1 = []
    y_deconv2 = []
    y_deconv3 = []
    y_raw1 = []
    y_raw2 = []
    y_raw3 = []
    inds = np.arange(n_responses)
    for n in range(1, n_responses+1):
        iters = n_responses // n
        np.random.shuffle(inds)
        for i in range(iters):
            tl = TraceList([deconv[k] for k in inds[n*i:n*(i+1)]])
            x.append(n)
            
            # calculate amplitude from averaged trace
            avg = tl.mean()
            y_deconv1.append(measure_amp(avg))
            
            # calculate average of amplitudes measured from individual traces
            amps = []
            for t in tl:
                amps.append(measure_amp(t))
            y_deconv2.append(np.mean(amps))
            
            # calculate a fake amplitude from the baseline region
            amps = []
            for t in tl:
                amps.append(measure_amp(t, baseline=(0, 2e-3), response=(6e-3, 10e-3)))
            y_deconv3.append(np.mean(amps))

            # calculate average amplitudes from raw traces
            tl = TraceList([responses.responses[k] for k in inds[n*i:n*(i+1)]])
            avg = tl.mean()
            y_raw1.append(measure_amp(avg))
            
            # calculate average of amplitudes measured from individual traces
            amps = []
            for t in tl:
                amps.append(measure_amp(t))
            y_raw2.append(np.mean(amps))

            # calculate a fake amplitude from the baseline region
            amps = []
            for t in tl:
                amps.append(measure_amp(t, baseline=(0, 2e-3), response=(6e-3, 10e-3)))
            y_raw3.append(np.mean(amps))
            

    ss = 7
    plt1 = pg.plot(labels={'bottom': 'n samples'}, title='Amplitudes from deconvolved traces')
    plt1.plot(x, y_deconv1, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 100), symbolSize=ss)
    plt1.plot(x, y_deconv2, pen=None, symbol='o', symbolPen=None, symbolBrush=(0, 255, 0, 100), symbolSize=ss)
    plt1.plot(x, y_deconv3, pen=None, symbol='x', symbolPen=None, symbolBrush=(255, 255, 0, 100), symbolSize=ss)
    
    plt2 = pg.plot(labels={'bottom': 'n samples'}, title='Amplitudes from raw traces')
    plt2.plot(x, y_raw1, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 0, 0, 100), symbolSize=ss)
    plt2.plot(x, y_raw2, pen=None, symbol='o', symbolPen=None, symbolBrush=(0, 255, 0, 100), symbolSize=ss)
    plt2.plot(x, y_raw3, pen=None, symbol='x', symbolPen=None, symbolBrush=(255, 255, 0, 100), symbolSize=ss)


    # Now see if the relationship between variance and n follows the predicted
    # 1/sqrt(n) curve
    for traces, plt in [(deconv, plt1), (responses.responses, plt2)]:
        all_amps = [measure_amp(t) for t in traces]
        stdev = np.std(all_amps)
        avg = np.mean(all_amps)
        nvals = np.arange(1, n_responses+1)
        pred = avg + stdev / nvals**0.5
        c1 = plt.plot(nvals, pred, pen='b')
        c2 = plt.plot(nvals, 2 * avg - pred, pen='b')
    
    
    
    