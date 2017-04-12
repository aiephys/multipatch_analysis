from copy import deepcopy
import numpy as np

from neuroanalysis.spike_detection import detect_evoked_spike
from neuroanalysis.stats import ragged_mean
from neuroanalysis.stimuli import square_pulses


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


class MultiPatchAnalyzer(Analyzer):
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
        
            [{pulse_n, pulse_ind, spike, response}, ...]
        
        """
        # detect presynaptic spikes
        pulse_stim = PulseStimAnalyzer.get(pre_rec)
        spikes = pulse_stim.evoked_spikes()
        
        # Select ranges to extract from postsynaptic recording
        result = []
        for i,pulse in enumerate(spikes):
            pulse = pulse.copy()
            spike = pulse['spike']
            if spike is None:
                continue
            
            # start recording window at the rising phase of the presynaptic spike
            pulse['rec_start'] = spike['rise_index']
            # stop 50 ms later
            pulse['rec_stop'] = spike['rise_index'] + int(50e-3 / pre_rec['primary'].dt)
            
            # truncate window early if there is another spike
            for pulse2 in spikes[i+1:]:
                if pulse2['spike'] is None:
                    continue
                pulse['rec_stop'] = min(pulse['rec_stop'], pulse2['spike']['rise_index'])
                break
            
            # Extract data from postsynaptic recording
            pulse['response'] = post_rec['primary'].data[pulse['rec_start']:pulse['rec_stop']]
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

    def get_evoked_responses(self, pre_id, post_id, clamp_mode='vc', stim_filter='20Hz'):
        all_spikes = self.all_evoked_responses()
        responses = []
        for rec in all_spikes[dev1][dev2]:
            
            # do filtering here:
            pre_rec = rec['pre_rec']
            post_rec = rec['post_rec']
            if post_rec.clamp_mode != 'vc':
                continue
            stim_name = pre_rec.meta['stim_name']
            if '20Hz' not in stim_name:
                continue
            
            
            for spike in rec['spikes']:
                if spike['spike'] is None:
                    continue
                responses.append(spike['response'])
        return responses
 
    def all_evoked_responses(self, max_freq=50):
        
        # loop over all sweeps (presynaptic)
        if self._all_spikes is None:
            all_spikes = {}
            for srec in self.expt.contents:
                mp_analyzer = MultiPatchAnalyzer(srec)
                
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


if __name__ == '__main__':
    #from experiment_list import ExperimentList
    #all_expts = ExperimentList(cache='expts_cache.pkl')

    ## pick one experiment with a lot of connections
    #for expt in all_expts:
        #if expt.expt_id[1] == '2017.03.20-0-0' and 'Pasha' in expt.expt_id[0]:
            #break

    import pyqtgraph as pg
    pg.dbg()

    from neuroanalysis.miesnwb import MiesNwb
    import sys
    expt_file = sys.argv[1]
    expt = MiesNwb(expt_file)
    
    analyzer = MultiPatchExperimentAnalyzer(expt)

    from neuroanalysis.ui.plot_grid import PlotGrid

    devs = analyzer.list_devs()
    n_devs = len(devs)
    plots = PlotGrid()
    plots.set_shape(n_devs, n_devs)
    plots.show() 
    
    for i, dev1 in enumerate(devs):
        for j, dev2 in enumerate(devs):
            if dev1 == dev2:
                continue
            responses = analyzer.get_evoked_responses(dev1, dev2)
            if len(responses) > 0:
                avg = ragged_mean(responses, method='clip')
                plots[i,j].plot(avg)
