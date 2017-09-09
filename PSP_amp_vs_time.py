"""
Question: how does PSP height change with sweep number.

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
import scipy.ndimage as ndi

from constants import INHIBITORY_CRE_TYPES
from constants import EXCITATORY_CRE_TYPES
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
    
    #expt_index = sys.argv[1]
    #pre_id, post_id = map(int, sys.argv[2:4])
    
    # Load experiment index
    cache_file = 'expts_cache.pkl'
    expts = ExperimentList(cache=cache_file)

    synapses = []
    for connection in expts.connection_summary():
        cells = connection['cells']
        expt = connection['expt']
        if cells[0].cre_type == 'pvalb' and cells[1].cre_type == 'pvalb':
            synapses.append((expt, cells[0].cell_id, cells[1].cell_id))
            
    #expt = expts[expt_index]
    
    #analyzer = MultiPatchExperimentAnalyzer(expt.data)
    #pulses = analyzer.get_evoked_responses(pre_id, post_id, clamp_mode='ic', pulse_ids=[0])
    p = pg.plot()
    
    for i,syn in enumerate(synapses):
        expt, pre_id, post_id = syn
        analyzer = DynamicsAnalyzer(expt, pre_id, post_id, align_to='spike')
        
        # collect all first pulse responses
        responses = analyzer.amp_group
        if len(responses) == 0:
            print("Skipping %s %d %d; no responses" % (expt.uid, pre_id, post_id))
            continue
        # collect all events
        #responses = analyzer.all_events
        
        n_responses = len(responses)
        
        # do exponential deconvolution on all responses
        deconv = TraceList()
#        grid1 = PlotGrid()
#        grid1.set_shape(2, 1)
        for j in range(n_responses):
            r = responses.responses[j]
#            grid1[0, 0].plot(r.time_values, r.data)
            
            filt = bessel_filter(r - np.median(r.time_slice(0, 10e-3).data), 300.)
            responses.responses[j] = filt
            
            dec = exp_deconvolve(r, 15e-3)
            baseline = np.median(dec.data[:100])
            r2 = bessel_filter(dec-baseline, 300.)
#            grid1[1, 0].plot(r2.time_values, r2.data)
            
            deconv.append(r2)
        
#        grid1.show()
        
    
        def measure_amp(trace, min_or_max, baseline=(6e-3, 8e-3), response=(13e-3, 17e-3)):
            baseline = trace.time_slice(*baseline).data.mean()
            if min_or_max=='max':
                peak = trace.time_slice(*response).data.max()
            elif min_or_max=='min':
                peak = trace.time_slice(*response).data.min()
            else:
                raise Exception('Are you looking for min or max')
            return peak - baseline
    
    
        average = responses.bsub_mean()
        max_peak = measure_amp(average, min_or_max='max', baseline=(6e-3, 8e-3), response=(12e-3, 16e-3))
        min_peak = measure_amp(average, min_or_max='min', baseline=(6e-3, 8e-3), response=(12e-3, 16e-3))
        
        max_min = "max" if abs(max_peak)> abs(min_peak) else "min"
        wrong_synapse_type_flag = False
        if max_min == "min" and expt.cells[pre_id].cre_type in EXCITATORY_CRE_TYPES: 
            print ("Whoa this synapse looks like a inhibitory when cre line would say it should be excitatory!!!" )  
            wrong_synapse_type_flag = True
        if max_min == "max" and expt.cells[pre_id].cre_type in INHIBITORY_CRE_TYPES: 
            print ("Whoa this synapse looks like a excitatory when cre line would say it should be inhibitory!!!" )    
            wrong_synapse_type_flag = True  
        peak=[]
        base=[]
        time=[]
        ordered=sorted(responses.responses, key=lambda rr:rr.start_time)
        for ii in range(len(responses)):
            rr=ordered[ii]
            peak.append(measure_amp(rr, min_or_max=max_min, baseline=(6e-3, 8e-3), response=(12e-3, 16e-3)))
            base.append(measure_amp(rr, min_or_max=max_min, baseline=(0e-3, 2e-3), response=(6e-3, 10e-3)))
            time.append(rr.start_time)
        mean_base=np.mean(base)
        time = np.array(time) - time[0]
        peak_minus_base_average=np.array(peak)-mean_base
        filtered=ndi.gaussian_filter(peak_minus_base_average, 2)    
        
#        p.plot(time, peak)
#        p.plot(time, base, pen='r')
#        p.plot(time, peak_minus_base_average, pen='g')
        if wrong_synapse_type_flag == False:
            p.plot(time, filtered, pen=(i, len(synapses)*1.3))
        elif wrong_synapse_type_flag == True:
            p.plot(time, filtered, pen='w')
        
        
            
        
        app.processEvents()
    
    