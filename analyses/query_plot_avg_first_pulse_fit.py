"""Queries and plots data from avg_first_pulse_fit table
"""

import numpy as np
import matplotlib.pyplot as plt
from multipatch_analysis.database import database as db
from multipatch_analysis import fit_average_first_pulse as fafp

#----------------------------------------------------------------
#--------------------------------query --------------------------
#----------------------------------------------------------------
'''Note that this query can take several minutes'''
session=db.Session()
data=session.query(fafp.AvgFirstPulseFit, db.Pair).join(db.Pair).all() #this need to correspond to import


for fit, pair in data:

    i_data=fit.ic_avg_psp_data
    v_data=fit.vc_avg_psp_data
    
    if (len(i_data) != 1) or len(v_data) != 1:
        plt.figure(figsize=(14,10))
        if len(i_data) != 1:  # DB contains an array with one value of 0 means there was no average wave to fit
            ax1=plt.subplot(1,1,1)
            scale_factor = 1e3 # converts to mV
            time = np.arange(len(i_data)) * fit.ic_dt * 1e3 # converting to ms
            ln1 = ax1.plot(time, i_data*scale_factor, 'b', label='current clamp')
            ln2 = ax1.plot(time, fit.ic_avg_psp_fit * scale_factor, 'r', label='nrmse=%f \namp (mV)=%f \nlatency (ms)=%f \nrise time (ms)=%f \ndecay tau=%f' % \
                                            (fit.ic_NRMSE, \
                                            fit.ic_amp * scale_factor, \
                                            fit.ic_latency*1e-3, \
                                            fit.ic_rise_time, \
                                            fit.ic_decay_tau))


            ax2=ax1.twinx()
            ln3=ax2.plot(time, fit.ic_weight, 'k', label='weight')
            ax1.set_ylabel('voltage (mv)')
            ax2.set_ylabel('weight')
            ax1.set_xlabel('time (ms): spike happens at 10 ms')

            lines_plot= ln1+ln2+ln3
            label_plot = [l.get_label() for l in lines_plot]
            ax1.legend(lines_plot, label_plot)
            plt.show()



    # elif clamp_mode == 'vc': 
    #     ln2=ax1.plot(waveform.time_values*1.e3, fit.best_fit*scale_factor, 'r', label='nrmse=%f \namp (pA)=%f \nlatency (ms)=%f \nrise time (ms)=%f \ndecay tau=%f' % \
    #                                     (fit.nrmse(), \
    #                                     fit.best_values['amp']*scale_factor, \
    #                                     (fit.best_values['xoffset']-time_before_spike)*1e3, \
    #                                     fit.best_values['rise_time']*1e3, \
    #                                     fit.best_values['decay_tau']))