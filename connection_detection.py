from copy import deepcopy
import numpy as np

from neuroanalysis.spike_detection import detect_evoked_spike



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
        if hasattr(rec, attr):
            raise TypeError("Object %s already has attached %s" % (rec, self.__class__.__name__))
        setattr(rec, attr, self)


class PulseStimAnalyzer(Analyzer):
    """Used for analyzing a patch clamp recording with square-pulse stimuli.
    """
    def __init__(self, rec):
        self._attach(rec)
        self.rec = rec
        self._pulse_inds = None
        self._evoked_spikes = None
        
    def pulse_inds(self):
        if self._pulse_inds is None:
            trace = self.rec['command'].data
            sdiff = np.diff(trace)
            on_inds = np.argwhere(sdiff > 0)[1:, 0]  # 1: skips test pulse
            off_inds = np.argwhere(sdiff < 0)[1:, 0]
            self._pulse_inds = on_inds, off_inds
        return self._pulse_inds

    def evoked_spikes(self):
        """Given presynaptic Recording, detect action potentials
        evoked by current injection or unclamped spikes evoked by a voltage pulse.
        """
        if self._evoked_spikes is None:
            pre_trace = self.rec['primary']

            # Detect pulse times
            pulses = zip(*self.pulse_inds())

            # detect spike times
            spike_info = []
            for i,pulse in enumerate(pulses):
                on, off = pulse
                spike = detect_evoked_spike(pre_rec, [on, off])
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

    def get_spike_responses(pre_rec, post_recs):
        """Given a presynaptic stimulus recording and a list of postsynaptic recordings,
        return a structure containing evoked responses from all postsynaptic recordings.
        
            [ {pre_rec, post_rec, [{pulse_n, pulse_ind, spike, response}, ...]}, ... ]
        
        """
        # detect presynaptic spikes
        spikes = detect_patch_evoked_spikes(pre_rec)
        
        # Select ranges to extract from postsynaptic recordings
        for i,pulse in enumerate(spikes):
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
        
        result = []
        
        
        raise Exception()
        #for rec in post_recs:
            #result.append(deepcopy(spikes))
            #for i,pulse in enumerate(result[-1]
            
        
        return result


def detect_connections(expt_data, max_freq=50):
    """
    loop over all sweeps (presynaptic)
        ignore sweeps with high induction frequency
        detect presynaptic spikes
        loop over all other sweeps in the same recording (postsynaptic)
            ignore sweeps that fail QC
            collect data following each stimulus
            collect baseline noise data
                (use noise.std() / sqrt(n) to estimate noise for averaged traces)
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
    
    # loop over all sweeps (presynaptic)
    for srec in expt_data.contents:
        print srec
        for pre_rec in srec.recordings:
            print "   ", pre_rec
            # todo: ignore sweeps with high induction frequency
            
            post_recs = srec.recordings[:]
            post_recs.remove(pre_rec)


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
    
    detect_connections(expt)
