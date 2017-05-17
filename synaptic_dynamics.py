import sys
import numpy as np
from connection_detection import MultiPatchSyncRecAnalyzer, EvokedResponseGroup, fit_psp
from database import ExperimentDatabase
from neuroanalysis.stats import ragged_mean
from neuroanalysis.baseline import float_mode
from neuroanalysis.ui.plot_grid import PlotGrid


if __name__ == '__main__':
    import pyqtgraph as pg
    from experiment_list import ExperimentList
    pg.mkQApp()
    pg.dbg()
    
    arg = sys.argv[1]
    expt_ind = int(arg)
    all_expts = ExperimentList(cache='expts_cache.pkl')
    expt = all_expts[expt_ind]

    pre = int(sys.argv[2])
    post = int(sys.argv[3])
    
    
    # Collect all data from NWB
    responses = {}
    for srec in expt.data.contents:
        if pre not in srec.devices or post not in srec.devices:
            continue
        pre_rec = srec[pre]
        post_rec = srec[post]
        if post_rec.clamp_mode != 'ic':
            continue
        
        analyzer = MultiPatchSyncRecAnalyzer.get(srec)
        resp = analyzer.get_spike_responses(pre_rec, post_rec)
        if len(resp) != 12:
            continue
        
        stim_params = analyzer.stim_params(pre_rec)
        responses.setdefault(stim_params, []).append(resp)

    # Sort responses by stimulus parameters and plot averages
    plots = PlotGrid()
    plots.set_shape(len(responses), 12)
    kinetics_group = EvokedResponseGroup()
    amp_group = EvokedResponseGroup()
    response_groups = {}
    for i,stim_params in enumerate(responses):
        resp = responses[stim_params]
        ind_freq, rec_delay = stim_params
        rec_delay = np.round(rec_delay, 2)
        response_groups[stim_params] = []
        for j in range(12):
            rg = EvokedResponseGroup()
            response_groups[stim_params].append(rg)
            for trial in resp:
                r = trial[j]['response']
                b = trial[j]['baseline']
                rg.add(r, b)
                plots[i,j].plot(r.time_values-r.time_values[0], r.data - np.median(b.data), pen=0.5)
                
                if ind_freq <= 50 and j in (7, 11):
                    kinetics_group.add(r, b)
                if ind_freq <= 100 and j == 0:
                    amp_group.add(r, b)
            avg = rg.bsub_mean()
            plots[i,j].plot(avg.time_values, avg.data, pen='g')
        plots[i,0].setLabels(left=("ind: %0.0f rec: %0.0f" % (ind_freq, rec_delay*1000), 'V'))
    plots.show()

    # Generate average first response
    avg_amp = amp_group.bsub_mean()
    amp_plot = pg.plot(title='First pulse amplitude')
    amp_plot.plot(avg_amp.time_values, avg_amp.data)

    # Make initial amplitude estimate
    ad = avg_amp.data
    dt = avg_amp.dt
    base = float_mode(ad[:int(10e-3/dt)])
    neg = ad[int(13e-3/dt):].min() - base
    pos = ad[int(13e-3/dt):].max() - base
    amp_est = neg if abs(neg) > abs(pos) else pos
    amp_plot.addLine(y=base + amp_est)
    amp_sign = '-' if amp_est < 0 else '+'

    # Generate average decay phase
    avg_kinetic = kinetics_group.bsub_mean()
    kin_plot = pg.plot(title='Kinetics')
    kin_plot.plot(avg_kinetic.time_values, avg_kinetic.data)
    
    # Make initial kinetics estimate
    kin_fit = fit_psp(avg_kinetic, sign=amp_sign, yoffset=0, amp=amp_est)
    kin_plot.plot(avg_kinetic.time_values, kin_fit.eval(), pen='b')
    
    # Fit all responses and plot dynamics curves
    with pg.ProgressDialog("Fitting responses..", maximum=len(response_groups)*12) as dlg:
        for i,stim_params in enumerate(response_groups):
            for j in range(12):
                rg = response_groups[stim_params][j]
                avg = rg.bsub_mean()
                fit = fit_psp(avg, sign=amp_sign, amp=amp_est, 
                            rise_time=kin_fit.best_values['rise_time'],
                            decay_tau=kin_fit.best_values['decay_tau'])
                plots[i,j].plot(avg.time_values, fit.eval(), pen='b')
                dlg += 1
                if dlg.wasCanceled():
                    raise Exception("Canceled response fit")
    





    
            






    #db = ExperimentDatabase()
    #db.load_data(expt, pre, post)

    #expts = db.get_table('experiment')
    #cell = db.get_table('cell')
    #srec = db.get_table('sync_rec')
    #rec = db.get_table('recording')
    #resp = db.get_table('response')
    #spike = db.get_table('stim_spike')
    #pulse = db.get_table('stim_pulse')

    ## don't need spikes for now; puse table gives us n_spikes
    ##spikes = pulse.merge(spike, how='left', left_on='id', right_on='pulse_id', suffixes=('_pulse', '_spike'))
    #pulses = pulse.merge(rec, how='left', left_on='recording_id', right_on='id', suffixes=('_pulse', '_pre_rec'))
    #pulses = pulses.merge(cell, how='left', left_on='cell_id', right_on='id', suffixes=('_spikes', '_pre_cell'))
    
    
    #resp = resp.merge(pulses, how='right', left_on='pulse_id', right_on='id_pulse', suffixes=('_resp', '_pulses'))

    #resp = rec.merge(resp, how='right', left_on='id', right_on='recording_id_resp', suffixes=('_post_rec', '_resp'))
    
    #resp = srec.merge(resp, how='right', left_on='id', right_on='sync_rec_id_post_rec', suffixes=('_srec', '_resp'))
    
    #resp = expts.merge(resp, how='right', left_on='id', right_on='expt_id_srec', suffixes=('_expt', '_resp'))
    
    #resp = resp.merge(cell, how='left', left_on='cell_id_post_rec', right_on='id', suffixes=('_resp', '_post_cell'))
    
    