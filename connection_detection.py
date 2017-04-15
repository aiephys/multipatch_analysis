from copy import deepcopy
import numpy as np
import pyqtgraph as pg

from neuroanalysis.spike_detection import detect_evoked_spike
from neuroanalysis.stats import ragged_mean
from neuroanalysis.stimuli import square_pulses
from neuroanalysis.data import Trace
from neuroanalysis.fitting import StackedPsp


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

    def get_spike_responses(self, pre_rec, post_rec):
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
        
        # select baseline region between 8th and 9th pulses
        # (ideally we should use more careful criteria for selecting a baseline region)
        dt = post_rec['primary'].dt
        baseline_dur = int(50e-3 / dt)
        stop = spikes[8]['pulse_ind'] - 1
        start = stop - baseline_dur
        baseline = post_rec['primary'][start:stop]
        
        # Select ranges to extract from postsynaptic recording
        result = []
        for i,pulse in enumerate(spikes):
            pulse = pulse.copy()
            spike = pulse['spike']
            if spike is None:
                continue
            
            # start recording window at the rising phase of the presynaptic spike
            pulse['rec_start'] = spike['rise_index'] - int(10e-3 / dt)
            
            if i+1 < len(spikes):
                # truncate window early if there is another pulse
                pulse['rec_stop'] = spikes[i+1]['pulse_ind']
            else:
                # otherwise, stop 50 ms later
                pulse['rec_stop'] = spike['rise_index'] + int(50e-3 / dt)
            
            # Extract data from postsynaptic recording
            pulse['response'] = post_rec['primary'][pulse['rec_start']:pulse['rec_stop']]
            
            # record pointer to the baseline before the first pulse
            pulse['baseline'] = baseline
            
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
        all_spikes = self.all_evoked_responses()
        responses = []
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
                    responses.append((resp, spike['baseline']))
        
        return responses
 
    def get_evoked_response_avg(self, pre_id, post_id, **kwds):
        """Return a baseline-subtracted, average evoked response trace between two cells.

        Source traces are collected using get_evoked_responses(), and all keyword arguments are
        forwarded there. All traces are downsampled to the minimum sample rate in the set.
        """
        responses = self.get_evoked_responses(pre_id, post_id, **kwds)
        if len(responses) == 0:
            return None

        # unzip into responses and corresponding baselines
        baselines = [r[1] for r in responses]
        responses = [r[0] for r in responses]
        
        # downsample all traces to the same rate
        # yarg: how does this change SNR?
        max_dt = max([trace.dt for trace in responses])
        downsampled = [trace.downsample(n=int(max_dt/trace.dt)).data for trace in responses]
        ds_baseline = [trace.downsample(n=int(max_dt/trace.dt)).data for trace in baselines]
        
        # average all together, clipping to minimum length
        avg = ragged_mean(downsampled, method='clip')
        avg_baseline = ragged_mean(ds_baseline, method='clip')
        
        # subtract baseline
        baseline = np.median(avg_baseline)
        bsub = avg - baseline

        result = Trace(bsub, dt=max_dt)

        # Attach some extra metadata to the result:
        result.meta['baseline'] = avg_baseline
        result.meta['baseline_med'] = baseline
        result.meta['baseline_std'] = avg_baseline.std()

        return result

    def all_evoked_responses(self):
        
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
        return self.all_evoked_responses().keys()


from neuroanalysis.ui.plot_grid import PlotGrid

def plot_response_averages(expt, **kwds):
    analyzer = MultiPatchExperimentAnalyzer.get(expt)
    devs = analyzer.list_devs()
    n_devs = len(devs)
    plots = PlotGrid()
    plots.set_shape(n_devs, n_devs)
    plots.show() 
    
    ranges = [([], []), ([], [])]
    points = []
    for i, dev1 in enumerate(devs):
        for j, dev2 in enumerate(devs):
            
            # select plot and adjust axes / labels
            plt = plots[i, j]
            plt.setXLink(plots[0, 0])
            plt.setYLink(plots[0, 0])
            plt.setLabels(bottom=(str(dev2), 's'))
            if kwds.get('clamp_mode', 'ic') == 'ic':
                plt.setLabels(left=('%s' % dev1, 'V'))
            else:
                plt.setLabels(left=('%s' % dev1, 'A'))

            if i < len(devs) - 1:
                plt.getAxis('bottom').setVisible(False)
            if j > 0:
                plt.getAxis('left').setVisible(False)            
            
            
            if dev1 == dev2:
                #plt.setVisible(False)
                continue
            
            print "==========", dev1, dev2
            avg_response = analyzer.get_evoked_response_avg(dev1, dev2, **kwds)
            if avg_response is not None:
                
                t = avg_response.time_values
                y = avg_response.data
                plt.plot(t, y, antialias=True)

                # fit!
                psp = StackedPsp()
                params = {
                    'xoffset': (10e-3, 9e-3, 20e-3),
                    'yoffset': (0, 'fixed'),
                    'rise_time': (10e-3, 0.1e-3, 30e-3),
                    'decay_tau': (50e-3, 1e-3, 200e-3),
                    'amp': (-1e-3, -100e-3, 100e-3),
                    'exp_amp': 'amp * amp_ratio',
                    'amp_ratio': (1, 0, 10),
                    'rise_power': (2, 'fixed'),
                }
                fit_kws = {'xtol': 1e-3, 'maxfev': 200, 'nan_policy': 'omit'}
                fit = psp.fit(y, x=t, fit_kws=fit_kws, params=params)
                print fit.best_values
                print "RMSE:", fit.rmse()
                nrmse = fit.nrmse()
                print "NRMSE:", nrmse
                snr = abs(fit.best_values['amp']) / avg_response.meta['baseline_std']
                print "SNR:", snr
                
                color = (
                    np.clip(255 * (1-nrmse), 0, 255),
                    np.clip(50 * snr, 0, 255),
                    np.clip(255 * nrmse, 0, 255)
                )

                plt.plot(t, fit.eval(), pen=color)


                points.append({'x': nrmse, 'y': snr, 'brush': color})
                # bl = avg_response.meta['baseline'] - avg_response.meta['baseline_med']
                # plt.plot(np.arange(len(bl)) * avg_response.dt, bl, pen='g')

                # keep track of data range across all plots
                ranges[0][0].append(y.min())
                ranges[0][1].append(y.max())
                ranges[1][0].append(t[0])
                ranges[1][1].append(t[-1])

    plots[0,0].setYRange(min(ranges[0][0]), max(ranges[0][1]))
    plots[0,0].setXRange(min(ranges[1][0]), max(ranges[1][1]))

    # scatter plot of SNR vs NRMSE
    plt = pg.plot()
    plt.setLabels(left='SNR', bottom='NRMSE')
    plt.plot([p['x'] for p in points], [p['y'] for p in points], pen=None, symbol='o', symbolBrush=[pg.mkBrush(p['brush']) for p in points])
        
    return plots



if __name__ == '__main__':
    #from experiment_list import ExperimentList
    #all_expts = ExperimentList(cache='expts_cache.pkl')

    ## pick one experiment with a lot of connections
    #for expt in all_expts:
        #if expt.expt_id[1] == '2017.03.20-0-0' and 'Pasha' in expt.expt_id[0]:
            #break

    import user
    import pyqtgraph as pg
    pg.dbg()

    from neuroanalysis.miesnwb import MiesNwb
    import sys
    expt_file = sys.argv[1]
    expt = MiesNwb(expt_file)
    
    plots = plot_response_averages(expt, clamp_mode='ic', min_duration=20e-3, pulse_ids=None)


"""
# run from connectivity_summary script:
import connection_detection as cd
cd.plot_response_averages(expts[171].data, min_duration=20e-3)

"""