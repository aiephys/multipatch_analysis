import sys
from collections import OrderedDict
import numpy as np
from connection_detection import MultiPatchSyncRecAnalyzer, EvokedResponseGroup, fit_psp
from neuroanalysis.stats import ragged_mean
from neuroanalysis.baseline import float_mode
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.fitting import PspTrain
from synaptic_release import ReleaseModel


def collect_stim_trains(expt, pre, post, padding):
    """Collect all stimulus-response recordings from an experiment between a
    specific pre- and post-synaptic cell.
    
    Returns data in 3 dicts, each keyed with the stimulus parameters (induction
    frequency and recovery delay):
        pulse_responses : {stim_params: [(pulse1, ..., pulse12), ...], ...}
            Postsynaptic responses separated into small chunks for each stimulus pulse
        train_responses : {stim_params: [(induction, recovery), ...], ...} 
            Postsynaptic responses, separated into larger chunks of multiple pulses for induction and recovery
        pulse_offsets : {stim_params: [pulse_offsets], ...}
            Offset times of pulses for each set of stimulus parameters
    """
    pre_pad, post_pad = padding
    pulse_responses = {}
    train_responses = {}
    pulse_offsets = {}
    for srec in expt.data.contents:
        if pre not in srec.devices or post not in srec.devices:
            continue
        pre_rec = srec[pre]
        post_rec = srec[post]
        if post_rec.clamp_mode != 'ic':
            continue
        
        analyzer = MultiPatchSyncRecAnalyzer.get(srec)
        resp = analyzer.get_spike_responses(pre_rec, post_rec, pre_pad=pre_pad)
        if len(resp) != 12:
            # for dynamics, we require all 12 pulses to elicit a presynaptic spike
            continue
        
        stim_params = analyzer.stim_params(pre_rec)
        pulse_responses.setdefault(stim_params, []).append(resp)
        
        ind, base = analyzer.get_train_response(pre_rec, post_rec, 0, 7, padding=(-pre_pad, post_pad))
        rec, base = analyzer.get_train_response(pre_rec, post_rec, 8, 11, padding=(-pre_pad, post_pad))
        ind.t0 = 0
        rec.t0 = 0
        if stim_params not in train_responses:
            train_responses[stim_params] = (EvokedResponseGroup(), EvokedResponseGroup())
        train_responses[stim_params][0].add(ind, base)
        train_responses[stim_params][1].add(rec, base)
        
        dt = pre_rec['command'].dt
        if stim_params not in pulse_offsets:
            i0 = resp[0]['pulse_ind']
            pulse_offsets[stim_params] = [(r['pulse_ind'] - i0)*dt for r in resp]
    
    # re-write as ordered dicts
    stim_param_order = sorted(pulse_offsets.keys())
    pulse_responses = OrderedDict([(k, pulse_responses[k]) for k in stim_param_order])
    train_responses = OrderedDict([(k, train_responses[k]) for k in stim_param_order])
    pulse_offsets = OrderedDict([(k, pulse_offsets[k]) for k in stim_param_order])
    
    return pulse_responses, train_responses, pulse_offsets


def get_kinetics_groups(pulse_responses):
    """Given a set of pulse responses, return EvokedResponseGroups that can be used to extract
    kinetic parameters.
    """
    kinetics_group = EvokedResponseGroup()
    amp_group = EvokedResponseGroup()
    for i,stim_params in enumerate(pulse_responses.keys()):
        # collect all individual pulse responses:
        #  - we can try fitting individual responses averaged across trials
        #  - collect first pulses for amplitude estimate
        #  - colect last pulses for kinetics estimate
        resp = pulse_responses[stim_params]
        ind_freq, rec_delay = stim_params
        rec_delay = np.round(rec_delay, 2)
        for j in range(12):
            for trial in resp:
                r = trial[j]['response']
                b = trial[j]['baseline']
                
                if ind_freq <= 20 or j in (7, 11):
                    kinetics_group.add(r, b)
                if ind_freq <= 100 and j == 0:
                    amp_group.add(r, b)
    
    return amp_group, kinetics_group


def plot_train_responses(train_responses):
    """
    Plot individual and averaged train responses for each set of stimulus parameters.
    
    Return a new PlotGrid.
    """
    train_plots = PlotGrid()
    train_plots.set_shape(len(stim_param_order), 2)

    for i,stim_params in enumerate(train_responses.keys()):
        
        # Collect and plot average traces covering the induction and recovery 
        # periods for this set of stim params
        ind_group = train_responses[stim_params][0]
        rec_group = train_responses[stim_params][1]
        
        for j in range(len(ind_group)):
            ind = ind_group.responses[j]
            rec = rec_group.responses[j]
            base = np.median(ind_group.baselines[j].data)
            train_plots[i,0].plot(ind.time_values, ind.data - base, pen=(128, 128, 128, 100))
            train_plots[i,1].plot(rec.time_values, rec.data - base, pen=(128, 128, 128, 100))
        ind_avg = ind_group.bsub_mean()
        rec_avg = rec_group.bsub_mean()

        ind_freq, rec_delay = stim_params
        train_plots[i,0].plot(ind_avg.time_values, ind_avg.data, pen='g', antialias=True)
        train_plots[i,1].plot(rec_avg.time_values, rec_avg.data, pen='g', antialias=True)
        train_plots[i,0].setLabels(left=("ind: %0.0f rec: %0.0f" % (ind_freq, rec_delay*1000), 'V'))
        
    train_plots.show()
    train_plots.setYLink(train_plots[0,0])
    for i in range(train_plots.shape[0]):
        train_plots[i,0].setXLink(train_plots[0,0])
        train_plots[i,1].setXLink(train_plots[0,1])
    train_plots.grid.ci.layout.setColumnStretchFactor(0, 3)
    train_plots.grid.ci.layout.setColumnStretchFactor(1, 2)
    train_plots.setClipToView(True)
    train_plots.setDownsampling(True, True, 'peak')
    
    return train_plots


if __name__ == '__main__':
    import pyqtgraph as pg
    from experiment_list import ExperimentList
    app = pg.mkQApp()
    pg.dbg()
    
    arg = sys.argv[1]
    expt_ind = int(arg)
    all_expts = ExperimentList(cache='expts_cache.pkl')
    expt = all_expts[expt_ind]

    # convert cell ID to headstage ID
    pre = int(sys.argv[2]) - 1
    post = int(sys.argv[3]) - 1

    pre_pad = 10e-3
    post_pad = 50e-3
    
    # Collect all data from NWB
    pulse_responses, train_responses, pulse_offsets = collect_stim_trains(expt, pre, post, padding=(pre_pad, post_pad))
    
    if len(pulse_responses) == 0:
        raise Exception("No suitable data found for cell %d -> cell %d in expt %s" % (pre, post, expt_ind))
           
    stim_param_order = pulse_offsets.keys()
    

    # Plot all individual and averaged train responses for all sets of stimulus parameters
    train_plots = plot_train_responses(train_responses)

    
    # Collect groups of events that can be averaged together to estimate the 
    # amplitude and kinetics of this synapse
    amp_group, kinetics_group = get_kinetics_groups(pulse_responses)
    
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
    avg_kinetic.t0 = 0
    kin_plot = pg.plot(title='Kinetics')
    kin_plot.plot(avg_kinetic.time_values, avg_kinetic.data)
    
    # Make initial kinetics estimate
    kin_fit = fit_psp(avg_kinetic, sign=amp_sign, yoffset=0, amp=amp_est, method='leastsq', fit_kws={})
    kin_plot.plot(avg_kinetic.time_values, kin_fit.eval(), pen='b')
    rise_time = kin_fit.best_values['rise_time']
    decay_tau = kin_fit.best_values['decay_tau']
    latency = kin_fit.best_values['xoffset'] - 10e-3

    ## Fit all responses and plot dynamics curves
    #with pg.ProgressDialog("Fitting responses..", maximum=len(response_groups)*12) as dlg:
        #for i,stim_params in enumerate(response_groups):
            #for j in range(12):
                #rg = response_groups[stim_params][j]
                #avg = rg.bsub_mean()
                #fit = fit_psp(avg, sign=amp_sign, amp=amp_est, 
                            #rise_time=(kin_fit.best_values['rise_time'], 'fixed'),
                            #decay_tau=(kin_fit.best_values['decay_tau'], 'fixed'))
                #plots[i,j].plot(avg.time_values, fit.eval(), pen='b')
                #dlg += 1
                #if dlg.wasCanceled():
                    #raise Exception("Canceled response fit")
    
    # Fit trains to multi-event models
    
    tasks = stim_param_order
    results = {}
    import pyqtgraph.multiprocess as mp
    with mp.Parallelize(enumerate(tasks), results=results, progressDialog='processing in parallel..') as tasker:
        for i,stim_params in tasker:
            grps = train_responses[stim_params]
            pulse_offset = pulse_offsets[stim_params]
            fits = []
            for j,grp in enumerate(grps):
                avg = grp.bsub_mean()
                
                base = np.median(avg.data[:int(10e-3/avg.dt)])
                
                # initial fit 
                
                args = {
                    'yoffset': (base, 'fixed'),
                    'xoffset': (0, -1e-3, 1e-3),
                    'rise_time': (rise_time, rise_time*0.5, rise_time*2),
                    'decay_tau': (decay_tau, decay_tau*0.5, decay_tau*2),
                    'rise_power': (2, 'fixed'),
                }
                
                pulses = [pulse_offset[:8], pulse_offset[8:]][j]
                for p,pt in enumerate(pulses):
                    args['xoffset%d'%p] = (pt - pulses[0] + pre_pad, 'fixed')
                    args['amp%d'%p] = (amp_est,) + tuple(sorted([0, amp_est * 10]))

                fit_kws = {'xtol': 1e-4, 'maxfev': 1000, 'nan_policy': 'omit'}                
                model = PspTrain(len(pulses))
                fit = model.fit(avg.data, x=avg.time_values, params=args, fit_kws=fit_kws, method='leastsq')
                

                # Fit again with decay tau per event
                # Slow, but might improve fit amplitudes
                args = {
                    'yoffset': (base, 'fixed'),
                    'xoffset': (0, -1e-3, 1e-3),
                    'rise_time': (fit.best_values['rise_time'], rise_time*0.5, rise_time*2),
                    'decay_tau': (fit.best_values['decay_tau'], decay_tau*0.5, decay_tau*2),
                    'rise_power': (2, 'fixed'),
                }
                
                for p,pt in enumerate(pulses):
                    args['xoffset%d'%p] = (fit.best_values['xoffset%d'%p], 'fixed')
                    args['amp%d'%p] = (fit.best_values['amp%d'%p],) + tuple(sorted([0, amp_est * 10]))
                    args['decay_tau_factor%d'%p] = (1, 0.5, 2)
                
                fit = model.fit(avg.data, x=avg.time_values, params=args, fit_kws=fit_kws, method='leastsq')
                
                fits.append((fit.best_values, len(pulses)))
                
            tasker.results[stim_params] = fits


    # plot train fits
    dyn_plots = PlotGrid()
    dyn_plots.set_shape(len(results), 1)
    spike_sets = []
    models = {4: PspTrain(4), 8: PspTrain(8)}
    for i,stim_params in enumerate(stim_param_order):
        fits = results[stim_params]
        amps = []
        for j,fit in enumerate(fits):
            fit, n_psp = fit
            print "-----------"
            print stim_params
            print fit
            import lmfit
            params = {k:lmfit.Parameter(name=k, value=v) for k,v in fit.items()}
            tvals = train_responses[stim_params][j].responses[0].time_values
            model = models[n_psp]
            train_plots[i,j].plot(tvals, model.eval(x=tvals, params=params), pen='b', antialias=True)
    
            amps.extend([abs(v) for k,v in sorted(fit.items()) if k.startswith('amp')])

        # plot dynamics
        dyn_plots[i,0].plot(amps)
        ind_freq, rec_delay = stim_params
        dyn_plots[i,0].setLabels(left=("ind: %0.0f rec: %0.0f" % (ind_freq, rec_delay*1000), 'V'))
            
        # prepare dynamics data for release model fit
        t = np.array(pulse_offsets[stim_params]) * 1000
        amps = np.array(amps) / amps[0]
        spike_sets.append((t, amps))
        
    dyn_plots.show()
    app.processEvents()
    
    # Fit release model to dynamics
    model = ReleaseModel()
    model.Dynamics['Dep'] = 1
    model.Dynamics['Fac'] = 1
    model.Dynamics['UR'] = 1
    model.Dynamics['SMR'] = 1
    model.Dynamics['DSR'] = 1

    fit = model.run_fit(spike_sets)
    rel_plots = PlotGrid()
    rel_plots.set_shape(2, 1)
    ind_plot = rel_plots[0, 0]
    ind_plot.setTitle('Release model fit - induction frequency')
    ind_plot.setLabels(bottom=('time', 's'), left='relative amplitude')
    rec_plot = rel_plots[1, 0]
    rec_plot.setTitle('Release model fit - recovery delay')
    rec_plot.setLabels(bottom=('time', 's'), left='relative amplitude')
    
    ind_plot.setLogMode(x=True, y=False)
    rec_plot.setLogMode(x=True, y=False)
    ind_plot.setXLink(rec_plot)
    
    for i,stim_params in enumerate(stim_param_order):
        x,y = spike_sets[i]
        output = model.eval(x, fit.values(), dt=0.5)
        y1 = output[:,1]
        x1 = output[:,0]
        if stim_params[1] == 0.252:
            ind_plot.plot((x+10)/1000., y, pen=None, symbol='o', symbolBrush=(i,10))
            ind_plot.plot((x1+10)/1000., y1, pen=(i,10))
        if stim_params[0] == 50:
            rec_plot.plot((x+10)/1000., y, pen=None, symbol='o', symbolBrush=(i,10))
            rec_plot.plot((x1+10)/1000., y1, pen=(i,10))
    
    rel_plots.show()






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
    
    