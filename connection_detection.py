from copy import deepcopy
import numpy as np
import scipy.signal
import pyqtgraph as pg

from neuroanalysis.spike_detection import detect_evoked_spike
from neuroanalysis.stats import ragged_mean
from neuroanalysis.stimuli import square_pulses
from neuroanalysis.data import Trace
from neuroanalysis.fitting import StackedPsp
from neuroanalysis.ui.plot_grid import PlotGrid


class Analyzer(object):
    @classmethod
    def get(cls, obj):
        """Get the analyzer attached to a recording, or create a new one.
        """
        analyzer = getattr(obj, '_' + cls.__name__, None)
        if analyzer is None:
            analyzer = cls(obj)
        return analyzer

    def _attach(self, obj):
        attr = '_' + self.__class__.__name__
        if hasattr(obj, attr):
            raise TypeError("Object %s already has attached %s" % (obj, self.__class__.__name__))
        setattr(obj, attr, self)


class PulseStimAnalyzer(Analyzer):
    """Used for analyzing a patch clamp recording with square-pulse stimuli.
    """
    def __init__(self, rec):
        self._attach(rec)
        self.rec = rec
        self._pulses = None
        self._evoked_spikes = None
        
    def pulses(self):
        """Return a list of (start, stop, amp) tuples describing square pulses
        in the stimulus.
        """
        if self._pulses is None:
            trace = self.rec['command'].data
            self._pulses = square_pulses(trace)
        return self._pulses

    def evoked_spikes(self):
        """Given presynaptic Recording, detect action potentials
        evoked by current injection or unclamped spikes evoked by a voltage pulse.
        """
        if self._evoked_spikes is None:
            pre_trace = self.rec['primary']

            # Detect pulse times
            pulses = self.pulses()

            # detect spike times
            spike_info = []
            for i,pulse in enumerate(pulses):
                on, off, amp = pulse
                if amp < 0:
                    # assume negative pulses do not evoke spikes
                    # (todo: should be watching for rebound spikes as well)
                    continue
                spike = detect_evoked_spike(self.rec, [on, off])
                spike_info.append({'pulse_n': i, 'pulse_ind': on, 'spike': spike})
            self._evoked_spikes = spike_info
        return self._evoked_spikes


class MultiPatchSyncRecAnalyzer(Analyzer):
    """Used for analyzing two or more synchronous patch clamp recordings where
    spikes are evoked in at least one and synaptic responses are recorded in
    others.
    """
    def __init__(self, srec):
        self._attach(srec)
        self.srec = srec

    def get_spike_responses(self, pre_rec, post_rec, align_to='pulse', pre_pad=10e-3):
        """Given a pre- and a postsynaptic recording, return a structure
        containing evoked responses.
        
            [{pulse_n, pulse_ind, spike, response, baseline}, ...]
        
        """
        # detect presynaptic spikes
        pulse_stim = PulseStimAnalyzer.get(pre_rec)
        spikes = pulse_stim.evoked_spikes()
        
        if len(spikes) < 10:
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
            if spike is None:
                continue
            
            if align_to == 'spike':
                # start recording window at the rising phase of the presynaptic spike
                pulse['rec_start'] = spike['rise_index'] - int(pre_pad / dt)
            elif align_to == 'pulse':
                # align to pulse onset
                pulse['rec_start'] = pulse['pulse_ind'] - int(pre_pad / dt)
            
            if i+1 < len(spikes):
                # truncate window early if there is another pulse
                pulse['rec_stop'] = spikes[i+1]['pulse_ind']
            else:
                # otherwise, stop 50 ms later
                pulse['rec_stop'] = spike['rise_index'] + int(50e-3 / dt)
            
            # Extract data from postsynaptic recording
            pulse['response'] = post_rec['primary'][pulse['rec_start']:pulse['rec_stop']]
            
            # select baseline region between 8th and 9th pulses
            # (ideally we should use more careful criteria for selecting a baseline region)
            baseline_dur = int(10e-3 / dt)
            stop = spikes[8]['pulse_ind'] - (i * baseline_dur)
            start = stop - baseline_dur
            pulse['baseline'] = post_rec['primary'][start:stop]
            
            result.append(pulse)
        
        return result


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
                    responses.add(resp, spike['baseline'])
        
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
            self.keys[i].add(key)

        data[item[-1]] = value
        self.keys[ndim-1].add(item[-1])

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
                
                t = avg_response.time_values
                y = avg_response.data
                plt.plot(t, y, antialias=True)

                # fit!                
                fit = responses[(dev1, dev2)].fit_psp(yoffset=0, mask_stim_artifact=(abs(dev1-dev2) < 3))
                lsnr = np.log(fit.snr)
                lerr = np.log(fit.nrmse())
                
                color = (
                    np.clip(255 * (-lerr/3.), 0, 255),
                    np.clip(50 * lsnr, 0, 255),
                    np.clip(255 * (1+lerr/3.), 0, 255)
                )

                plt.plot(t, fit.best_fit, pen=color)
                # plt.plot(t, fit.init_fit, pen='y')

                points.append({'x': lerr, 'y': lsnr, 'brush': color})

                if show_baseline:
                    # plot baseline for reference
                    bl = avg_response.meta['baseline'] - avg_response.meta['baseline_med']
                    plt.plot(np.arange(len(bl)) * avg_response.dt, bl, pen=(0, 100, 0), antialias=True)

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
    def __init__(self, pre_id, post_id, **kwds):
        self.pre_id = pre_id
        self.post_id = post_id
        self.kwds = kwds
        self.responses = []
        self.baselines = []
        self._bsub_mean = None

    def add(self, response, baseline):
        self.responses.append(response)
        self.baselines.append(baseline)
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
            avg = trace_mean(responses)
            avg_baseline = trace_mean(baselines).data
            max_dt = avg.dt

            # subtract baseline
            baseline = np.median(avg_baseline)
            bsub = avg.data - baseline

            result = Trace(bsub, dt=max_dt)

            # Attach some extra metadata to the result:
            result.meta['baseline'] = avg_baseline
            result.meta['baseline_med'] = baseline
            result.meta['baseline_std'] = scipy.signal.detrend(avg_baseline).std()

            self._bsub_mean = result

        return self._bsub_mean

    def fit_psp(self, **kwds):
        response = self.bsub_mean()
        if response is None:
            return None
        return fit_psp(response, **kwds)


def trace_mean(traces):
    """Return the mean of a list of traces.

    Downsamples to the minimum rate and clips ragged edges.
    """
    max_dt = max([trace.dt for trace in traces])
    downsampled = [trace.downsample(n=int(max_dt/trace.dt)).data for trace in traces]
    avg = ragged_mean(downsampled, method='clip')
    return Trace(avg, dt=max_dt)





def fit_psp(response, mode='ic', sign='any', xoffset=11e-3, yoffset=(0, 'fixed'), mask_stim_artifact=True):
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
    
    params = []
    for amp, amp_min, amp_max in amps:
        p2 = base_params.copy()
        p2['amp'] = (amp, amp_min, amp_max)
        params.append(p2)

    dt = response.dt
    weight = np.ones(len(y))
    weight[:int(10e-3/dt)] = 0.5
    if mask_stim_artifact:
        # Use zero weight for fit region around the stimulus artifact
        weight[int(10e-3/dt):int(13e-3/dt)] = 0

    fit_kws = {'weights': weight, 'xtol': 1e-4, 'maxfev': 300, 'nan_policy': 'omit'}
    
    best_fit = None
    best_score = None
    for p in params:
        fit = psp.fit(y, x=t, params=p, fit_kws=fit_kws, method='leastsq')
        err = np.sum(fit.residual**2)
        if best_fit is None or err < best_score:
            best_fit = fit
            best_score = err
    fit = best_fit

    # nrmse = fit.nrmse()
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


if __name__ == '__main__':

    import user
    import pyqtgraph as pg
    pg.dbg()

    from neuroanalysis.miesnwb import MiesNwb
    import sys
    arg = sys.argv[1]
    try:
        from experiment_list import ExperimentList
        all_expts = ExperimentList(cache='expts_cache.pkl')
        expt_ind = int(arg)
        expt = all_expts[expt_ind].data
    except ValueError:
        expt_file = arg
        expt = MiesNwb(expt_file)
    
    plots = plot_response_averages(expt, show_baseline=True, clamp_mode='ic', min_duration=15e-3, pulse_ids=None)

    detect_connections(expt)
