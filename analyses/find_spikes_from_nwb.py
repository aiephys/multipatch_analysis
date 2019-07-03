"""prototyping spike detection on file from nwb.
"""
import multipatch_analysis.database as db
import multipatch_analysis.connection_strength as cs
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from weighted_error_lib import * 
from neuroanalysis.data import Trace, TraceList
import pdb
#import ipfx.feature_extractor as fe
import ipfx.spike_detector as sd
from scipy.optimize import curve_fit
from multipatch_analysis.data import Analyzer, PulseStimAnalyzer



#import pyqtgraph as pg #this is here to be able to use pyqt debugger 
#pg.dbg() #will open console if exception happens and you are running in interactive mode


# cells with currently poorly identified spikes 
cell_ids = [#[1544582617.589, 1, 8, 5656957, 5654136, 5654045],  #this is a good text bc two fail but the others are sort of sad looking.
            #[1544582617.589, 1, 6, 5654136, 5654045], 
#            [1497417667.378, 5, 2, 7483977, 7483912],
#            [1491942526.646, 8, 1, 6693052, 6693000],  #this one is giving me issues #this presynaptic cell is sick.  Spiking is ambiguious, very interesting examples
            [1521004040.059, 5, 6],
            [1534293227.896, 7, 8, 7271530], #this one is not quite perfected from the recurve up at the end fix with derivat
            [1540356446.981, 8, 6],
#            [1550101654.271, 1, 6], # these spike toward the end and are found correctly
            [1516233523.013, 6, 7],  #very interesting example: a voltage deflection happens very early but cant be seen in dvvdt due to to onset being to early.  Think about if there is a way to fix this.  Maybe and initial pulse window.  
            [1534297702.068, 7, 2]
            ]

s = db.Session()
for cell_id in cell_ids:
    expt = db.experiment_from_timestamp(cell_id[0])
    pair = expt.pairs[cell_id[1], cell_id[2]]
    synapse = pair.synapse
    synapse_type = pair.connection_strength.synapse_type
    ic_pulse_ids = pair.avg_first_pulse_fit.ic_pulse_ids
    pulse_responses = pair.pulse_responses

    # get unique sweeps from pulse responses
    sweeps_in_pr = []
    for pr in pulse_responses:
#        if pr.stim_pulse_id in ic_pulse_ids:# cell_id[3:]: # 
        if pr.stim_pulse.recording.patch_clamp_recording.clamp_mode == 'ic':
            sweep = pr.stim_pulse.recording.sync_rec.ext_id
            channel = pr.stim_pulse.recording.electrode.device_id
            sweeps_in_pr.append((sweep, channel)) 
    unique_sweeps = list(set(sweeps_in_pr))

    # run spike detector on sweeps grabbed from nwb
    for sweep, channel in unique_sweeps: 
        pre_rec = expt.data.contents[sweep][channel]
        pulse_stim = PulseStimAnalyzer.get(pre_rec)
        spikes = pulse_stim.evoked_spikes()

        # view on a plot
        plt.figure(figsize =  (10, 8))
        voltage_mv = pre_rec['primary'].data
        dvdt = np.diff(voltage_mv) 
        d2vdt2 = np.diff(dvdt)
        time_ms = pre_rec['primary'].time_values
        ax1=plt.subplot(111)
        ax1.plot(time_ms, voltage_mv)
        ax2=ax1.twinx()
        ax2.plot(time_ms[1:], dvdt, color='r')
        ax2.plot(time_ms[2:], d2vdt2, color='g')
        for s in spikes:
            if s['spike']:
                plt.axvline(s['spike']['max_dvdt_time'], color='k')
        plt.show()
    
