import sys
from collections import OrderedDict
import numpy as np
from connection_detection import MultiPatchSyncRecAnalyzer, EvokedResponseGroup, fit_psp
from neuroanalysis.stats import ragged_mean
from neuroanalysis.baseline import float_mode
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.fitting import PspTrain
from neuroanalysis.synaptic_release import ReleaseModel


class DynamicsAnalyzer(object):
    def __init__(self, expt, pre_cell, post_cell):
        self.expt = expt
        self.pre_cell = pre_cell
        self.post_cell = post_cell

        # how much padding to include when extracting events
        self.pre_pad = 10e-3
        self.post_pad = 50e-3
        
        self._reset()
        
    def _reset(self):
        """Clear out cached analysis results
        """
        self._pulse_responses = None
        self._train_responses = None
        self._pulse_offsets = None
        
        self._amp_group = None
        self._kinetics_group = None
        
        self._psp_estimate = {}
        
        self._train_fit_results = None
        self._spike_sets = None
        
        self._last_model_fit = None

    @property
    def pulse_responses(self):
        if self._pulse_responses is None:
            self._collect_stim_trains()
        return self._pulse_responses

    @property
    def train_responses(self):
        if self._train_responses is None:
            self._collect_stim_trains()
        return self._train_responses

    @property
    def pulse_offsets(self):
        if self._pulse_offsets is None:
            self._collect_stim_trains()
        return self._pulse_offsets
    
    @property
    def stim_param_order(self):
        return self.pulse_responses.keys()

    @property
    def amp_group(self):
        if self._amp_group is None:
            self._get_kinetics_groups()
        return self._amp_group

    @property
    def kinetics_group(self):
        if self._kinetics_group is None:
            self._get_kinetics_groups()
        return self._kinetics_group

    @property
    def psp_estimate(self):
        if 'amp' not in self._psp_estimate:
            self.estimate_amplitude()
            
        if 'rise_time' not in self._psp_estimate:
            self.estimate_kinetics()
            
        return self._psp_estimate

    @property
    def train_fit_results(self):
        if self._train_fit_results is None:
            self.fit_response_trains()
        return self._train_fit_results
    
    @property
    def spike_sets(self):
        if self._spike_sets is None:
            self.prepare_spike_sets()
        return self._spike_sets
    
    def _collect_stim_trains(self):
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
        expt = self.expt
        
        # convert cell ID to headstage ID
        pre = self.pre_cell-1
        post = self.post_cell-1
        
        pre_pad, post_pad = self.pre_pad, self.post_pad
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
        
        self._pulse_responses = pulse_responses
        self._train_responses = train_responses
        self._pulse_offsets = pulse_offsets

    def _get_kinetics_groups(self):
        """Given a set of pulse responses, return EvokedResponseGroups that can be used to extract
        kinetic parameters.
        """
        pulse_responses = self.pulse_responses
        kinetics_group = EvokedResponseGroup()
        amp_group = EvokedResponseGroup()
        for i,stim_params in enumerate(pulse_responses.keys()):
            # collect all individual pulse responses:
            #  - we can try fitting individual responses averaged across trials
            #  - collect first pulses for amplitude estimate
            #  - colect last pulses for kinetics estimate
            resp = pulse_responses[stim_params]
            ind_freq, rec_delay = stim_params
            for j in range(12):
                for trial in resp:
                    r = trial[j]['response']
                    b = trial[j]['baseline']
                    
                    if ind_freq <= 20 or j in (7, 11):
                        kinetics_group.add(r, b)
                    if ind_freq <= 100 and j == 0:
                        amp_group.add(r, b)
        
        self._amp_group = amp_group
        self._kinetics_group = kinetics_group

    def plot_train_responses(self, plot_grid=None):
        """
        Plot individual and averaged train responses for each set of stimulus parameters.
        
        Return a new PlotGrid.
        """
        train_responses = self.train_responses
        
        if plot_grid is None:
            train_plots = PlotGrid()
        else:
            train_plots = plot_grid
        train_plots.set_shape(len(self.stim_param_order), 2)

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
            rec_delay = np.round(rec_delay, 2)
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

    def estimate_amplitude(self, plot=False):
        amp_group = self.amp_group
        amp_est = None
        amp_plot = None
        amp_sign = None
        avg_amp = None
        if len(amp_group) == 0:
            return amp_est, amp_sign, avg_amp, amp_plot
        # Generate average first response
        avg_amp = amp_group.bsub_mean()
        if plot:
            amp_plot = pg.plot(title='First pulse amplitude')
            amp_plot.plot(avg_amp.time_values, avg_amp.data)

        # Make initial amplitude estimate
        ad = avg_amp.data
        dt = avg_amp.dt
        base = float_mode(ad[:int(10e-3/dt)])
        neg = ad[int(13e-3/dt):].min() - base
        pos = ad[int(13e-3/dt):].max() - base
        amp_est = neg if abs(neg) > abs(pos) else pos
        if plot:
            amp_plot.addLine(y=base + amp_est)
        amp_sign = '-' if amp_est < 0 else '+'
        
        self._psp_estimate['amp'] = amp_est
        self._psp_estimate['amp_sign'] = amp_sign
        
        return amp_est, amp_sign, avg_amp, amp_plot

    def estimate_kinetics(self, plot=False):
        kinetics_group = self.kinetics_group
        
        # Generate average decay phase
        avg_kinetic = kinetics_group.bsub_mean()
        avg_kinetic.t0 = 0
        
        if plot:
            kin_plot = pg.plot(title='Kinetics')
            kin_plot.plot(avg_kinetic.time_values, avg_kinetic.data)
        else:
            kin_plot = None
        
        # Make initial kinetics estimate
        kin_fit = fit_psp(avg_kinetic, sign=amp_sign, yoffset=0, amp=amp_est, method='leastsq', fit_kws={})
        if plot:
            kin_plot.plot(avg_kinetic.time_values, kin_fit.eval(), pen='b')
        rise_time = kin_fit.best_values['rise_time']
        decay_tau = kin_fit.best_values['decay_tau']
        latency = kin_fit.best_values['xoffset'] - 10e-3

        self._psp_estimate['rise_time'] = rise_time
        self._psp_estimate['decay_tau'] = decay_tau
        self._psp_estimate['latency'] = latency
        
        return rise_time, decay_tau, latency, kin_plot

    def fit_response_trains(self):
        train_responses = self.train_responses
        pulse_offsets = self._pulse_offsets
        
        tasks = train_responses.keys()
        results = OrderedDict([(task,None) for task in tasks])
        import pyqtgraph.multiprocess as mp
        with mp.Parallelize(enumerate(tasks), results=results, progressDialog='Fitting PSP trains..') as tasker:
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
                        args['xoffset%d'%p] = (pt - pulses[0] + self.pre_pad, 'fixed')
                        args['amp%d'%p] = (amp_est,) + tuple(sorted([0, amp_est * 10]))

                    fit_kws = {'xtol': 1e-4, 'maxfev': 3000, 'nan_policy': 'omit'}                
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

        self._train_fit_results = results
        return results

    def prepare_spike_sets(self):
        """Generate spike amplitude structure needed for release model fitting
        """
        train_responses = self.train_responses
        pulse_offsets = self.pulse_offsets
        results = self.train_fit_results
        
        spike_sets = []
        for i,stim_params in enumerate(results.keys()):
            fits = results[stim_params]
            amps = []
            for j,fit in enumerate(fits):
                fit, n_psp = fit
                amps.extend([abs(v) for k,v in sorted(fit.items()) if k.startswith('amp')])

            # prepare dynamics data for release model fit
            t = np.array(pulse_offsets[stim_params]) * 1000
            amps = np.array(amps) / amps[0]
            spike_sets.append((t, amps))
            
        self._spike_sets = spike_sets

    def plot_train_fits(self, plot_grid):
        train_responses = self.train_responses
        pulse_offsets = self.pulse_offsets
        results = self.train_fit_results
        train_plots = plot_grid

        #dyn_plots = PlotGrid()
        #dyn_plots.set_shape(len(results), 1)
        models = {4: PspTrain(4), 8: PspTrain(8)}
        for i,stim_params in enumerate(results.keys()):
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
        
                #amps.extend([abs(v) for k,v in sorted(fit.items()) if k.startswith('amp')])

            # plot dynamics
            #dyn_plots[i,0].plot(amps)
            #ind_freq, rec_delay = stim_params
            #dyn_plots[i,0].setLabels(left=("ind: %0.0f rec: %0.0f" % (ind_freq, rec_delay*1000), 'V'))
            
        #dyn_plots.show()

    def fit_release_model(self, dynamics):
        model = ReleaseModel()
        for gate in dynamics:
            if gate not in model.Dynamics:
                raise ValueError("Unknown gating mechanism for release model: %s" % gate)
        for gate in model.Dynamics:
            if gate in dynamics:
                model.Dynamics[gate] = 1
                
        fit = model.run_fit(self.spike_sets)
        
        self._last_model_fit = (model, fit)
        return model, fit

    def plot_model_results(self, model=None, fit=None):
        if model is None:
            if self._last_model_fit is None:
                raise Exception("Must run fit_release_model before plotting results.")
            model, fit = self._last_model_fit
        spike_sets = self.spike_sets
        
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
        
        for i,stim_params in enumerate(self.stim_param_order):
            x,y = spike_sets[i]
            output = model.eval(x, fit.values(), dt=0.5)
            y1 = output[:,1]
            x1 = output[:,0]
            if stim_params[1] - 0.250 < 5e-3:
                ind_plot.plot((x+10)/1000., y, pen=None, symbol='o', symbolBrush=(i,10))
                ind_plot.plot((x1+10)/1000., y1, pen=(i,10))
            if stim_params[0] == 50:
                rec_plot.plot((x+10)/1000., y, pen=None, symbol='o', symbolBrush=(i,10))
                rec_plot.plot((x1+10)/1000., y1, pen=(i,10))
        
        rel_plots.show()
        return rel_plots



if __name__ == '__main__':
    import pyqtgraph as pg
    from experiment_list import ExperimentList
    app = pg.mkQApp()
    pg.dbg()
    
    arg = sys.argv[1]
    expt_ind = int(arg)
    all_expts = ExperimentList(cache='expts_cache.pkl')
    expt = all_expts[expt_ind]

    pre_cell = int(sys.argv[2])
    post_cell = int(sys.argv[3])


    analyzer = DynamicsAnalyzer(expt, pre_cell, post_cell)
    if len(analyzer.pulse_responses) == 0:
        raise Exception("No suitable data found for cell %d -> cell %d in expt %s" % (pre_cell, post_cell, expt_ind))
           
    # Plot all individual and averaged train responses for all sets of stimulus parameters
    train_plots = analyzer.plot_train_responses()

    if '--no-fit' in sys.argv:
        sys.exit(0)  # user requested no fitting; bail out early

    # Estimate PSP amplitude
    amp_est, amp_sign, amp_plot = analyzer.estimate_amplitude(plot=True)
    app.processEvents()
    
    # Estimate PSP kinetics
    rise_time, decay_tau, latency, kin_plot = analyzer.estimate_kinetics(plot=True)
    app.processEvents()
    
    # Fit trains to multi-event models and plot the results
    analyzer.plot_train_fits(train_plots)

    # update GUI before doing model fit
    app.processEvents()

    if '--no-model' in sys.argv:
        sys.exit(0)  # user requested no model; bail out early
    
    # Fit release model to dynamics
    model, fit = analyzer.fit_release_model(dynamics=['Dep', 'Fac', 'UR'])
    
    # plot nost recently-generated fit results
    rel_plots = analyzer.plot_model_results()






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
    
    