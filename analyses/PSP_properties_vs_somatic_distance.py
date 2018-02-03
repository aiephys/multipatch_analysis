"""
Question:  How do synaptic properties vary with somatic distance 
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
import matplotlib.pyplot as plt

import datetime
from multipatch_analysis.constants import INHIBITORY_CRE_TYPES
from multipatch_analysis.constants import EXCITATORY_CRE_TYPES
from multipatch_analysis.connection_detection import MultiPatchExperimentAnalyzer, EvokedResponseGroup
from multipatch_analysis.synaptic_dynamics import DynamicsAnalyzer
from multipatch_analysis.experiment_list import cached_experiments
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import TraceList, PatchClampRecording
from neuroanalysis.filter import bessel_filter
from neuroanalysis.event_detection import exp_deconvolve
from scipy import stats
import allensdk.core.json_utilities as ju

def check_synapse(expt, cells):
    '''checks if a synapse meets the requirements and if so, it appends it to the synapse 
    dictionary.  A synapses list must be initialized before calling this function
    inputs:
        expt: object
            object obtained from cached_experiments.connection_summary[*]['expt']
        cells: object
            object obtained from cached_experiments.connection_summary[*]['cells']    
    output:
        returns nothing but appends info to the synapse list
    '''
    try: #needed because pyqt is breaking on datetime sometimes
        if expt.expt_info['solution']=='2mM Ca & Mg':
            synapses.append((expt, cells[0].cell_id, cells[1].cell_id))
    except:
        pass


if __name__ == '__main__':
    app = pg.mkQApp()
    pg.dbg()
    
    # Load experiment index
    expts = cached_experiments()
#    expts.select(calcium='high')  #this is throwing datetime errors

    connection_list=[
                     ['rorb', 'rorb'],
                     ['tlx3', 'tlx3'],
                     ['ntsr1', 'ntsr1'],
                     ['L23pyr', 'L23pyr'],
                     ['sim1','sim1']]

    dictionary={}
    for synapic_pairs in connection_list:
        print(synapic_pairs)
    
        synapses = []
        for connection in expts.connection_summary():
            cells = connection['cells']
            expt = connection['expt']
            pre_synaptic=synapic_pairs[0]
            post_synaptic=synapic_pairs[1]
            if pre_synaptic=='L23pyr':
                if cells[0].target_layer=='2/3' and cells[1].target_layer=='2/3':
                    check_synapse(expt, cells)
            else:
                if cells[0].cre_type == pre_synaptic and cells[1].cre_type == post_synaptic:
                    check_synapse(expt, cells)
        
#        title_str= pre_synaptic+' to '+post_synaptic
#        time_vs_psp_plot = pg.plot(labels={'left': 'peak of synaptic deflection (V)', 
#                            'bottom': 'time since first recorded synapse (s)', 
#                            'top':(title_str+' connections: progression of synaptic defection over an experiment')})    
#        ave_psp_plot = pg.plot(labels={'top':('average base-line subtracted first pulse synaptic deflection ('+ title_str+ ')'), 
#                          'bottom': 'time (s)', 
#                          'left':'voltage (V)'}) 
#        sweep_vs_psp_plot = pg.plot(labels={'left': 'peak of synaptic deflection (V)', 
#                            'bottom': 'sweep number', 
#                            'top':(title_str+' connections: progression of synaptic defection over an experiment')})  
#        
        raw=[]
        filtered=[]
        time_list=[]
        sweep_number_list=[]
        num_of_synapses=0
        for i,syn in enumerate(synapses):
            expt, pre_id, post_id = syn
            analyzer = DynamicsAnalyzer(expt, pre_id, post_id, align_to='spike')
            
            # collect all first pulse responses
            amp_responses = analyzer.amp_group
            if len(amp_responses) == 0:
                print("Skipping %s %d %d; no responses" % (expt.uid, pre_id, post_id))
                continue
            
            EvokedResponseGroup(amp_responses)
            ERG=EvokedResponseGroup(pre_id, post_id).amp_responses
            model=ERG.fit_psp()
            print (model)