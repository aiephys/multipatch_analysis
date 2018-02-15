import numpy as np

from neuroanalysis.miesnwb import MiesNwb, MiesSyncRecording, MiesRecording
from neuroanalysis.stimuli import square_pulses
from neuroanalysis.spike_detection import detect_evoked_spike


class MultiPatchExperiment(MiesNwb):
    """Extension of neuroanalysis data abstraction layer to include
    multipatch-specific metadata.
    """
    def create_sync_recording(self, sweep_id):
        return MultiPatchSyncRecording(self, sweep_id)

        
class MultiPatchSyncRecording(MiesSyncRecording):
    def __init__(self, nwb, sweep_id):
        MiesSyncRecording.__init__(self, nwb, sweep_id)
        self._baseline_mask = None
        try:
            self.meta['temperature'] = self.recordings[0].meta['notebook']['Async AD 1: Bath Temperature']
        except Exception:
            self.meta['temperature'] = None
    
    def create_recording(self, sweep_id, ch):
        miesrec = MiesRecording(self, sweep_id, ch)
        stim = miesrec.meta['stim_name'].lower()
        if 'pulsetrain' in stim or 'recovery' in stim:
            return MultiPatchProbe(miesrec)
        else:
            return miesrec

    def baseline_regions(self, settling_time=100e-3):
        """Return a list of start,stop pairs indicating regions during the recording that are expected to be quiescent
        due to absence of pulses.
        """
        if self._baseline_mask is None:
            mask = np.zeros(len(self.recordings[0]['primary']), dtype=bool)
            dt = self.recordings[0]['primary'].dt
            settle_size = int(settling_time / dt)
            for rec in self.recordings:
                pa = PulseStimAnalyzer.get(rec)
                for pulse in pa.pulses():
                    mask[pulse[0]:pulse[1] + settle_size] = True
            self._baseline_mask = mask

            starts = list(np.argwhere(~mask[1:] & mask[:-1])[:,0])
            stops = list(np.argwhere(mask[1:] & ~mask[:-1])[:,0])
            if starts[0] > stops[0]:
                starts.insert(0, 0)
            if stops[-1] < starts[-1]:
                stops.append(len(mask))
            self._baseline_regions = zip(starts, stops)

        return self._baseline_regions


class MultiPatchProbe(MiesRecording):
    def __init__(self, recording):
        self._parent_rec = recording
        self._base_regions = None
        
    #@property    
    #def induction_frequency(self):
        #return self.stim_params
    def __len__(self):
        return len(self._parent_rec)

    def __getattr__(self, attr):
        if '_parent_rec' not in self.__dict__:
            raise AttributeError(attr)
        return getattr(self._parent_rec, attr)

    @property
    def baseline_regions(self):
        # detect baseline regions from pulses rather than labnotebook
        # (can't always count on there being delay periods in the recording)
        if self._base_regions is None:
            self._base_regions = self._parent_rec.parent.baseline_regions()
        return self._base_regions


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

    def stim_params(self):
        """Return induction frequency and recovery delay.
        """
        pulses = [p[0] for p in self.pulses() if p[2] > 0]
        if len(pulses) < 2:
            return None, None
        dt = self.rec['command'].dt
        ind_freq = np.round(1.0 / (dt * (pulses[1] - pulses[0])))
        rec_delay = np.round(dt*np.diff(pulses).max(), 3)
        
        return ind_freq, rec_delay
