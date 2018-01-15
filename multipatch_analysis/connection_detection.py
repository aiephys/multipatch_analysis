from copy import deepcopy
import numpy as np
import scipy.signal
import pyqtgraph as pg

from .data import MultiPatchProbe, Analyzer, PulseStimAnalyzer
from neuroanalysis.stats import ragged_mean
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.fitting import StackedPsp
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.filter import bessel_filter




class BaselineDistributor(Analyzer):
    """Used to find baseline regions in a trace and distribute them on request.
    """
    def __init__(self, rec):
        self._attach(rec)
        self.rec = rec
        self.baselines = rec.baseline_regions
        self.ptr = 0

    def get_baseline_chunk(self, duration=20e-3):
        while True:
            if len(self.baselines) == 0:
                return None
            bl = self.baselines[0]
            remain = bl.dt * (len(bl) - self.ptr)
            if remain >= duration:
                size = int(duration / bl.dt)
                chunk = bl[self.ptr:self.ptr+size]
                self.ptr += size
                return chunk
            else:
                self.baselines.pop(0)


class MultiPatchSyncRecAnalyzer(Analyzer):
    """Used for analyzing two or more synchronous patch clamp recordings where
    spikes are evoked in at least one and synaptic responses are recorded in
    others.
    """
    def __init__(self, srec):
        self._attach(srec)
        self.srec = srec

    def get_spike_responses(self, pre_rec, post_rec, align_to='pulse', pre_pad=10e-3, require_spike=True):
        """Given a pre- and a postsynaptic recording, return a structure
        containing evoked responses.
        
            [{pulse_n, pulse_ind, spike, response, baseline}, ...]
        
        """
        # detect presynaptic spikes
        pulse_stim = PulseStimAnalyzer.get(pre_rec)
        spikes = pulse_stim.evoked_spikes()
        
        baseline_dist = BaselineDistributor.get(post_rec)

        if not isinstance(pre_rec, MultiPatchProbe):
            # this does not look like the correct kind of stimulus; bail out
            # Ideally we can make this agnostic to the exact stim type in the future,
            # but for now we rely on the delay period between pulses 8 and 9 to get
            # a baseline measurement.
            return []

        dt = post_rec['primary'].dt
        
        # Select ranges to extract from postsynaptic recording
        result = []
        for i,pulse in enumerate(spikes):
            pulse = pulse.copy()
            spike = pulse['spike']
            if require_spike and spike is None:
                continue
            
            if align_to == 'spike':
                # start recording window at the rising phase of the presynaptic spike
                pulse['rec_start'] = spike['rise_index'] - int(pre_pad / dt)
            elif align_to == 'pulse':
                # align to pulse onset
                pulse['rec_start'] = pulse['pulse_ind'] - int(pre_pad / dt)
            
            max_stop = pulse['rec_start'] + int(50e-3 / dt)
            if i+1 < len(spikes):
                # truncate window early if there is another pulse
                pulse['rec_stop'] = min(max_stop, spikes[i+1]['pulse_ind'])
            else:
                # otherwise, stop 50 ms later
                pulse['rec_stop'] = max_stop
            
            # Extract data from postsynaptic recording
            pulse['response'] = post_rec['primary'][pulse['rec_start']:pulse['rec_stop']]

            # Extract presynaptic spike and stimulus command
            pulse['pre_rec'] = pre_rec['primary'][pulse['rec_start']:pulse['rec_stop']]
            pulse['command'] = pre_rec['command'][pulse['rec_start']:pulse['rec_stop']]

            # select baseline region between 8th and 9th pulses
            baseline_dur = int(100e-3 / dt)
            stop = spikes[8]['pulse_ind']
            start = stop - baseline_dur
            pulse['baseline'] = post_rec['primary'][start:stop]
            pulse['baseline_start'] = start
            pulse['baseline_stop'] = stop
            assert len(pulse['baseline']) > 0

            result.append(pulse)
        
        return result

    def get_pulse_response(self, pre_rec, post_rec, first_pulse=0, last_pulse=-1):
        pulse_stim = PulseStimAnalyzer.get(pre_rec)
        spikes = pulse_stim.evoked_spikes()
        
        if len(spikes) < 10:
            return None
        
        dt = post_rec['primary'].dt
        start = spikes[first_pulse]['pulse_ind'] - int(20e-3 / dt)
        stop = spikes[last_pulse]['pulse_ind'] + int(50e-3 / dt)
        return post_rec['primary'][start:stop]

    def stim_params(self, pre_rec):
        return PulseStimAnalyzer.get(pre_rec).stim_params()
        
    def get_train_response(self, pre_rec, post_rec, start_pulse, stop_pulse, padding=(-10e-3, 50e-3)):
        """Return the part of the post-synaptic recording during a range of pulses,
        along with a baseline chunk
        """
        pulse_stim = PulseStimAnalyzer.get(pre_rec)
        pulses = [p[0] for p in pulse_stim.pulses() if p[2] > 0]
        
        post_trace = post_rec['primary']
        pre_trace = pre_rec['primary']
        stim_command = pre_rec['command']

        dt = post_trace.dt
        start = pulses[start_pulse] + int(padding[0]/dt)
        stop = pulses[stop_pulse] + int(padding[1]/dt)
        assert start > 0
        
        response = post_trace[start:stop]
        pre_spike = pre_trace[start:stop]
        command = stim_command[start:stop]
        
        bstop = pulses[8] - int(10e-3/dt)
        bstart = bstop - int(200e-3/dt)
        baseline = post_trace[bstart:bstop]
        
        return response, baseline, pre_spike, command

    def find_artifacts(self, rec, freq=None, threshold=-10e-3):
        """Find contaminating artifacts in recording such as spontaneous spikes or electrical noise, return True
        state if anything is found"""

        artifacts = False
        # find anything that crosses a max threshold of -10mV
        if np.any(rec >= threshold):
            artifacts = True
        return artifacts
        


class MultiPatchExperimentAnalyzer(Analyzer):
    """
    loop over all sweeps (presynaptic)
        ignore sweeps with high induction frequency
        detect presynaptic spikes
        loop over all other sweeps in the same recording (postsynaptic)
            ignore sweeps that fail QC
            collect data following each stimulus
            collect baseline noise data
                use noise.std() / sqrt(n) to estimate noise for averaged traces
                use mode(noise) to correct for baseline
            collect mode / holding
    Average together all collected sweep chunks from the same cell, and of the same holding and clamp mode
        all pulses
        pulses 1, 2, 9, 10
        pulses 7, 8, 11, 12
    Fit result to psp
    Compare fit amplitude to noise floor
    Z-score fit NRMSE against fits from pre-stimulus data
        correct for multiple comparisons?
    Additional metric to detect low release probability connections?
    """
    def __init__(self, expt):
        self._attach(expt)
        self.expt = expt
        self._all_spikes = None

    def get_evoked_responses(self, pre_id, post_id, clamp_mode='ic', stim_filter=None, min_duration=None, pulse_ids=None):
        """Return all evoked responses from device pre_id to post_id with the given
        clamp mode and stimulus conditions.
        
        All traces are *downsampled* to the minimum sample rate across the set
        of returned responses.
        
        Returns a list of (response, baseline) pairs. 
        """
        all_spikes = self._all_evoked_responses()
        responses = EvokedResponseGroup(pre_id, post_id)
        for rec in all_spikes[pre_id][post_id]:
            
            # do filtering here:
            pre_rec = rec['pre_rec']
            post_rec = rec['post_rec']
            
            if post_rec.clamp_mode != clamp_mode:
                continue
            
            if stim_filter is not None:
                stim_name = pre_rec.meta['stim_name']
                if stim_filter not in stim_name:
                    continue
            
            for spike in rec['spikes']:
                if spike['spike'] is None:
                    continue
                if pulse_ids is not None and spike['pulse_n'] not in pulse_ids:
                    continue
                resp = spike['response']
                if resp.duration >= min_duration:
                    responses.add(resp.copy(t0=0), spike['baseline'])
        
        return responses
 
    def get_evoked_response_matrix(self, **kwds):
        """Returned evoked responses for all pre/post pairs
        """
        devs = self.list_devs()

        all_responses = {}
        rows = set()
        cols = set()
        for i, dev1 in enumerate(devs):
            for j, dev2 in enumerate(devs):
                if dev1 == dev2:
                    continue
                resp = self.get_evoked_responses(dev1, dev2, **kwds)
                if len(resp) > 0:
                    rows.add(dev1)
                    cols.add(dev2)
                all_responses[(dev1, dev2)] = resp
        rows = sorted(list(rows))
        cols = sorted(list(cols))

        return all_responses, rows, cols

    def _all_evoked_responses(self):
        
        # loop over all sweeps (presynaptic)
        if self._all_spikes is None:
            all_spikes = {}
            for srec in self.expt.contents:
                mp_analyzer = MultiPatchSyncRecAnalyzer.get(srec)
                
                for pre_rec in srec.recordings:
                    if not isinstance(pre_rec, MultiPatchProbe):
                        continue
                    pre_id = pre_rec.device_id
                    all_spikes.setdefault(pre_id, {})
                    # todo: ignore sweeps with high induction frequency
                    
                    for post_rec in srec.recordings:
                        if post_rec is pre_rec:
                            continue
                        post_id = post_rec.device_id
                        all_spikes[pre_id].setdefault(post_id, [])
                        spikes = {
                            'spikes': mp_analyzer.get_spike_responses(pre_rec, post_rec),
                            'pre_rec': pre_rec,
                            'post_rec': post_rec,
                        }
                        all_spikes[pre_id][post_id].append(spikes)
            self._all_spikes = all_spikes
        return self._all_spikes

    def list_devs(self):
        return self._all_evoked_responses().keys()



class NDDict(object):
    """N-dimensional dictionary.
    """
    def __init__(self, ndim):
        self._ndim = int(ndim)
        self._data = {}
        self._keys = [set() for i in range(ndim)]

    @property
    def ndim(self):
        return self._ndim

    def __setitem__(self, item, value):
        ndim = self.ndim
        assert len(item) == ndim
        data = self._data
        for i,key in enumerate(item[:-1]):
            data = data.setdefault(key, {})
            self.keys()[i].add(key)

        data[item[-1]] = value
        self.keys()[ndim-1].add(item[-1])

    def keys(self):
        return self._keys



def plot_response_averages(expt, show_baseline=False, **kwds):
    analyzer = MultiPatchExperimentAnalyzer.get(expt)
    devs = analyzer.list_devs()

    # First get average evoked responses for all pre/post pairs
    responses, rows, cols = analyzer.get_evoked_response_matrix(**kwds)

    # resize plot grid accordingly
    plots = PlotGrid()
    plots.set_shape(len(rows), len(cols))
    plots.show() 
    
    ranges = [([], []), ([], [])]
    points = []

    # Plot each matrix element with PSP fit
    for i, dev1 in enumerate(rows):
        for j, dev2 in enumerate(cols):
            # select plot and hide axes
            plt = plots[i, j]
            if i < len(devs) - 1:
                plt.getAxis('bottom').setVisible(False)
            if j > 0:
                plt.getAxis('left').setVisible(False)

            if dev1 == dev2:
                plt.getAxis('bottom').setVisible(False)
                plt.getAxis('left').setVisible(False)
                continue
            
            # adjust axes / labels
            plt.setXLink(plots[0, 0])
            plt.setYLink(plots[0, 0])
            plt.addLine(x=10e-3, pen=0.3)
            plt.addLine(y=0, pen=0.3)
            plt.setLabels(bottom=(str(dev2), 's'))
            if kwds.get('clamp_mode', 'ic') == 'ic':
                plt.setLabels(left=('%s' % dev1, 'V'))
            else:
                plt.setLabels(left=('%s' % dev1, 'A'))

            
            # print "==========", dev1, dev2
            avg_response = responses[(dev1, dev2)].bsub_mean()
            if avg_response is not None:
                avg_response.t0 = 0
                t = avg_response.time_values
                y = bessel_filter(Trace(avg_response.data, dt=avg_response.dt), 2e3).data
                plt.plot(t, y, antialias=True)

                # fit!                
                #fit = responses[(dev1, dev2)].fit_psp(yoffset=0, mask_stim_artifact=(abs(dev1-dev2) < 3))
                #lsnr = np.log(fit.snr)
                #lerr = np.log(fit.nrmse())
                
                #color = (
                    #np.clip(255 * (-lerr/3.), 0, 255),
                    #np.clip(50 * lsnr, 0, 255),
                    #np.clip(255 * (1+lerr/3.), 0, 255)
                #)

                #plt.plot(t, fit.best_fit, pen=color)
                ## plt.plot(t, fit.init_fit, pen='y')

                #points.append({'x': lerr, 'y': lsnr, 'brush': color})

                #if show_baseline:
                    ## plot baseline for reference
                    #bl = avg_response.meta['baseline'] - avg_response.meta['baseline_med']
                    #plt.plot(np.arange(len(bl)) * avg_response.dt, bl, pen=(0, 100, 0), antialias=True)

                # keep track of data range across all plots
                ranges[0][0].append(y.min())
                ranges[0][1].append(y.max())
                ranges[1][0].append(t[0])
                ranges[1][1].append(t[-1])

    plots[0,0].setYRange(min(ranges[0][0]), max(ranges[0][1]))
    plots[0,0].setXRange(min(ranges[1][0]), max(ranges[1][1]))

    # scatter plot of SNR vs NRMSE
    plt = pg.plot()
    plt.setLabels(left='ln(SNR)', bottom='ln(NRMSE)')
    plt.plot([p['x'] for p in points], [p['y'] for p in points], pen=None, symbol='o', symbolBrush=[pg.mkBrush(p['brush']) for p in points])
    # show threshold line
    line = pg.InfiniteLine(pos=[0, 6], angle=180/np.pi * np.arctan(1))
    plt.addItem(line, ignoreBounds=True)

    return plots


class EvokedResponseGroup(object):
    """A group of similar synaptic responses.

    This is intended to be used as a container for many repeated responses evoked from
    a single pre/postsynaptic pair. It provides methods for computing the average,
    baseline-subtracted response and for fitting the average to a curve.
    """
    def __init__(self, pre_id=None, post_id=None, **kwds):
        self.pre_id = pre_id
        self.post_id = post_id
        self.kwds = kwds
        self.responses = []
        self.baselines = []
        self.spikes = []
        self.commands = []
        self._bsub_mean = None

    def add(self, response, baseline, pre_spike=None, stim_command=None):
        self.responses.append(response)
        self.baselines.append(baseline)
        self.spikes.append(pre_spike)
        self.commands.append(stim_command)
        self._bsub_mean = None

    def __len__(self):
        return len(self.responses)

    def bsub_mean(self):
        """Return a baseline-subtracted, average evoked response trace between two cells.

        All traces are downsampled to the minimum sample rate in the set.
        """
        if len(self) == 0:
            return None

        if self._bsub_mean is None:
            responses = self.responses
            baselines = self.baselines
            
            # downsample all traces to the same rate
            # yarg: how does this change SNR?

            avg = TraceList([r.copy(t0=0) for r in responses]).mean()
            avg_baseline = TraceList([b.copy(t0=0) for b in baselines]).mean().data

            # subtract baseline
            baseline = np.median(avg_baseline)
            bsub = avg.data - baseline
            result = avg.copy(data=bsub)
            assert len(result.time_values) == len(result)

            # Attach some extra metadata to the result:
            result.meta['baseline'] = avg_baseline
            result.meta['baseline_med'] = baseline
            if len(avg_baseline) == 0:
                result.meta['baseline_std'] = None
            else:
                result.meta['baseline_std'] = scipy.signal.detrend(avg_baseline).std()

            self._bsub_mean = result

        return self._bsub_mean

    def mean(self):
        if len(self) == 0:
            return None
        return TraceList(self.responses).mean()

    def fit_psp(self, **kwds):
        response = self.bsub_mean()
        if response is None:
            return None
        return fit_psp(response, **kwds)


def fit_psp(response, mode='ic', sign='any', xoffset=11e-3, yoffset=(0, 'fixed'), mask_stim_artifact=True, method='leastsq', fit_kws=None, **kwds):
    t = response.time_values
    y = response.data

    if mode == 'ic':
        amp = .2e-3
        amp_max = 100e-3
        rise_time = 5e-3
        decay_tau = 50e-3
    elif mode == 'vc':
        amp = 20e-12
        amp_max = 500e-12
        rise_time = 1e-3
        decay_tau = 4e-3
    else:
        raise ValueError('mode must be "ic" or "vc"')

    amps = [(amp, 0, amp_max), (-amp, -amp_max, 0)]
    if sign == '-':
        amps = amps[1:]
    elif sign == '+':
        amps = amps[:1]
    elif sign != 'any':
        raise ValueError('sign must be "+", "-", or "any"')

    psp = StackedPsp()
    base_params = {
        'xoffset': (xoffset, 10e-3, 15e-3),
        'yoffset': yoffset,
        'rise_time': (rise_time, rise_time/2., rise_time*2.),
        'decay_tau': (decay_tau, decay_tau/10., decay_tau*10.),
        'exp_amp': 'amp * amp_ratio',
        'amp_ratio': (1, 0, 10),
        'rise_power': (2, 'fixed'),
    }
    
    if 'rise_time' in kwds:
        rt = kwds.pop('rise_time')
        if not isinstance(rt, tuple):
            rt = (rt, rt/2., rt*2.)
        base_params['rise_time'] = rt
                
    if 'decay_tau' in kwds:
        dt = kwds.pop('decay_tau')
        if not isinstance(dt, tuple):
            dt = (dt, dt/2., dt*2.)
        base_params['decay_tau'] = dt
                
    base_params.update(kwds)
    
    params = []
    for amp, amp_min, amp_max in amps:
        p2 = base_params.copy()
        p2['amp'] = (amp, amp_min, amp_max)
        params.append(p2)

    dt = response.dt
    weight = np.ones(len(y))
    #weight[:int(10e-3/dt)] = 0.5
    if mask_stim_artifact:
        # Use zero weight for fit region around the stimulus artifact
        weight[int(10e-3/dt):int(12e-3/dt)] = 0
    weight[int(12e-3/dt):int(19e-3/dt)] = 30

    if fit_kws is None:
        fit_kws = {'xtol': 1e-4, 'maxfev': 300, 'nan_policy': 'omit'}

    if 'weight' not in fit_kws:
        fit_kws['weights'] = weight
    
    best_fit = None
    best_score = None
    for p in params:
        fit = psp.fit(y, x=t, params=p, fit_kws=fit_kws, method=method)
        err = np.sum(fit.residual**2)
        if best_fit is None or err < best_score:
            best_fit = fit
            best_score = err
    fit = best_fit

    # nrmse = fit.nrmse()
    if 'baseline_std' in response.meta:
        fit.snr = abs(fit.best_values['amp']) / response.meta['baseline_std']
        fit.err = fit.rmse() / response.meta['baseline_std']
    # print fit.best_values
    # print "RMSE:", fit.rmse()
    # print "NRMSE:", nrmse
    # print "SNR:", snr

    return fit


def detect_connections(expt):
    analyzer = MultiPatchExperimentAnalyzer.get(expt)

    # First get average evoked responses for all pre/post pairs with long decay time
    all_responses, rows, cols = analyzer.get_evoked_response_matrix(clamp_mode='ic', min_duration=16e-3)

    for pre_id in rows:
        for post_id in cols:
            try:
                response = all_responses[(pre_id, post_id)]
                if len(response) == 0:
                    continue
            except KeyError:
                continue

            # fit average to extract PSP decay
            fit = response.fit_psp(yoffset=0)

            # make connectivity call
            lsnr = np.log(fit.snr)
            lnrmse = np.log(fit.nrmse())
            if lsnr > lnrmse + 6:
                print "Connection:", pre_id, post_id, fit.snr, fit.nrmse()

