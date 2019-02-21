"""Queries and plots voltage and current clamp data fits from avg_first_pulse_fit table.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from multipatch_analysis.database import database as db
from multipatch_analysis import fit_average_first_pulse as fafp

#----------------------------------------------------------------
#-------------------------- queries -----------------------------
#----------------------------------------------------------------
session=db.Session()

# ---- Query all---------------
'''Note that this query can take several minutes'''
#data=session.query(fafp.AvgFirstPulseFit, db.Pair).join(db.Pair).all() #this need to correspond to import

# --- Query some good PV examples--------
pv = [(1533244490.755, 6, 4),
    (1530559621.966, 7, 6),
    (1527020350.517, 6, 7),
    (1525903868.115, 6, 5),
    (1525903868.115, 7, 6),
    (1525812130.134, 7, 6),
    (1523398687.484, 2, 1),
    (1523398687.484, 2, 3),
    (1523398687.484, 3, 2),
    (1521667891.153, 3, 4),
    (1519425770.224, 6, 7),
    (1519425770.224, 7, 2),
    (1518564082.242, 4, 6),
    (1517356361.727, 7, 6),
    (1517356063.393, 3, 4),
    (1517348193.989, 1, 6),
    (1517348193.989, 6, 8),
    (1517266234.424, 5, 4),
    (1517266234.424, 6, 5),
    (1517269966.858, 2, 4),
    (1517269966.858, 2, 8),
    (1517269966.858, 7, 8),
    (1516820219.855, 3, 4),
    (1516744107.347, 3, 6),
    (1513976168.029, 2, 1),
    (1511913615.871, 5, 4),
    (1510268339.831, 4, 5),
    (1510268339.831, 4, 7),
    (1510268339.831, 5, 4),
    (1510268339.831, 7, 4),
    (1508189280.042, 4, 8),
    (1507235159.776, 5, 6),
    (1505421643.243, 6, 7),
    (1505421643.243, 7, 6),
    (1505515553.146, 3, 8),
    (1505515553.146, 6, 5),
    (1505515553.146, 6, 7),
    (1505515553.146, 7, 3),
    (1500668871.652, 1, 7),
    (1496703451.204, 4, 1),
    (1495574744.012, 4, 2),
    (1493159586.902, 3, 5),
    (1493159586.902, 7, 8),
    (1492018873.073, 8, 4),
    (1491587831.024, 8, 3),
    (1491252584.113, 1, 5),
    (1491347329.805, 3, 6),
    (1491347329.805, 6, 2),
    (1490043541.218, 2, 7),
    (1490043541.218, 6, 2),
    (1490043541.218, 7, 2),
    (1490043541.218, 7, 6),
    (1484952266.115, 8, 4),
    (1539801888.714, 1, 8),
    (1539801888.714, 7, 1)]

data = []
for uid, pre_cell_ext_id, post_cell_ext_id in pv:

    pre_cell = db.aliased(db.Cell)
    post_cell = db.aliased(db.Cell)
    stuff = session.query(fafp.AvgFirstPulseFit, db.Pair)\
            .join(db.Pair).join(db.Experiment)\
            .join(pre_cell, db.Pair.pre_cell_id==pre_cell.id)\
            .join(post_cell, db.Pair.post_cell_id==post_cell.id)\
            .filter(db.Experiment.acq_timestamp == uid)\
            .filter(pre_cell.ext_id == pre_cell_ext_id)\
            .filter(post_cell.ext_id == post_cell_ext_id).all()  
    data.append([stuff[0][0], stuff[0][1]])  # this format is used by plotting below


#------------------------------------------------------------------------------
#-----------plot voltage and current clamp from query output ------------------
#------------------------------------------------------------------------------
for fit, pair in data:
    """data: list of list containing and pair.  
            i.e. if there is only one entry it should look like [[fit, pair]]
    """

    i_data=fit.ic_avg_psp_data
    v_data=fit.vc_avg_psp_data
    
    title='%0.3f, cells %i %s to %i %s; distance=%.1f um' % (pair.experiment.acq_timestamp, pair.pre_cell.ext_id, pair.pre_cell.cre_type, pair.post_cell.ext_id, pair.post_cell.cre_type, pair.distance*1e6)
    if (len(i_data) != 1) or len(v_data) != 1:
        plt.figure(figsize=(14,10))
        ax1=plt.subplot(2,1,1)
        ax3=plt.subplot(2,1,2)

        if len(i_data) != 1:  # DB contains an array with one value of 0 means there was no average wave to fit
            # plotting current clamp
            scale_factor = 1e3 # converts to mV
            time = np.arange(len(i_data)) * fit.ic_dt * 1e3 # converting to ms
            ln1 = ax1.plot(time, i_data*scale_factor, 'b', label='current clamp, n=%i' % len(fit.ic_pulse_ids))
            ln2 = ax1.plot(time, fit.ic_avg_psp_fit * scale_factor, 'r', label='nrmse=%f \namp (mV)=%f \nlatency (ms)=%f \nrise time (ms)=%f \ndecay tau=%f' % \
                                            (fit.ic_nrmse, \
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

        # plot voltage clamp
        if len(v_data) != 1:  # DB contains an array with one value of 0 means there was no average wave to fit
            scale_factor = 1.e12 # converts to pA
            time = np.arange(len(v_data)) * fit.vc_dt * 1e3 # converting to ms
            ln4 = ax3.plot(time, v_data*scale_factor, 'b', label='voltage clamp, n=%i' % len(fit.vc_pulse_ids))
            ln5 = ax3.plot(time, fit.vc_avg_psp_fit * scale_factor, 'r', label='nrmse=%f \namp (pA)=%f \nlatency (ms)=%f \nrise time (ms)=%f \ndecay tau=%f' % \
                                            (fit.vc_nrmse, \
                                            fit.vc_amp * scale_factor, \
                                            fit.vc_latency*1e-3, \
                                            fit.vc_rise_time, \
                                            fit.vc_decay_tau))


            ax4=ax3.twinx()
            ln6=ax4.plot(time, fit.vc_weight, 'k', label='weight')
            ax3.set_ylabel('current (pA)')
            ax4.set_ylabel('weight')
            ax3.set_xlabel('time (ms): spike happens at 10 ms')

            lines_plot= ln4+ln5+ln6
            label_plot = [l.get_label() for l in lines_plot]
            ax3.legend(lines_plot, label_plot)

        ax1.set_title(title, fontsize = 20)    
        plt.show()
