import numpy as np

from neuroanalysis.miesnwb import MiesNwb, MiesSyncRecording, MiesRecording
from neuroanalysis.stimuli import find_square_pulses
from neuroanalysis.spike_detection import detect_evoked_spikes


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
        stim = miesrec.meta['notebook']['Stim Wave Name'].lower()
        if 'pulsetrain' in stim or 'recovery' in stim:
            return MultiPatchProbe(miesrec)
        else:
            return miesrec

    def baseline_regions(self, settling_time=100e-3):
        """Return a list of start,stop pairs indicating regions during the recording that are expected to be quiescent
        due to absence of pulses.
        """
        if self._baseline_mask is None:
            pri = self.recordings[0]['primary']
            mask = np.zeros(len(pri), dtype=bool)
            dt = pri.dt
            settle_size = int(settling_time / dt)
            for rec in self.recordings:
                pa = PulseStimAnalyzer.get(rec)
                for pulse in pa.pulses():
                    start = pri.index_at(pulse[0])
                    stop = pri.index_at(pulse[1])
                    mask[start:stop + settle_size] = True
            self._baseline_mask = mask

            starts = list(np.argwhere(~mask[1:] & mask[:-1])[:,0])
            stops = list(np.argwhere(mask[1:] & ~mask[:-1])[:,0])
            if starts[0] > stops[0]:
                starts.insert(0, 0)
            if stops[-1] < starts[-1]:
                stops.append(len(mask))
            baseline_inds = [r for r in zip(starts, stops) if r[1] > r[0]]
            self._baseline_regions = [(pri.time_at(i0), pri.time_at(i1)) for i0, i1 in baseline_inds]

        return self._baseline_regions


class MultiPatchProbe(MiesRecording):
    def __init__(self, recording):
        self._parent_rec = recording
        self._base_regions = None
        
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
        """Return a list of (start_time, stop_time, amp) tuples describing square pulses
        in the stimulus.
        """
        if self._pulses is None:
            trace = self.rec['command']
            pulses = find_square_pulses(trace)
            self._pulses = []
            for p in pulses:
                start = p.global_start_time
                stop = p.global_start_time + p.duration
                self._pulses.append((start, stop, p.amplitude))
        return self._pulses

    def pulse_chunks(self):
        """Return time-slices of this recording where evoked spikes are expected to be found (one chunk
        per pulse)
        
        Each recording returned has extra metadata keys added: 
        - pulse_edges: start/end times of the stimulus pulse
        - pulse_amplitude: amplitude of stimulus puse (in V or A)
        - pulse_n: the number of this pulse (all detected square pulses are numbered in order from 0)

        """
        pre_trace = self.rec['primary']

        # Detect pulse times
        pulses = self.pulses()

        # cut out a chunk for each pulse
        chunks = []
        for i,pulse in enumerate(pulses):
            pulse_start_time, pulse_end_time, amp = pulse
            if amp < 0:
                # assume negative pulses do not evoke spikes
                # (todo: should be watching for rebound spikes as well)
                continue
            # cut out a chunk of the recording for spike detection
            start_time = pulse_start_time - 2e-3
            stop_time = pulse_end_time + 4e-3
            if i < len(pulses) - 1:
                # truncate chunk if another pulse is present
                next_pulse_time = pulses[i+1][0]
                stop_time = min(stop_time, next_pulse_time)
            chunk = self.rec.time_slice(start_time, stop_time)
            chunk.meta['pulse_edges'] = [pulse_start_time, pulse_end_time]
            chunk.meta['pulse_amplitude'] = amp
            chunk.meta['pulse_n'] = i
            chunks.append(chunk)
        return chunks

    def evoked_spikes(self):
        """Given presynaptic Recording, detect action potentials
        evoked by current injection or unclamped spikes evoked by a voltage pulse.
        """
        if self._evoked_spikes is None:
            spike_info = []
            for i,chunk in enumerate(self.pulse_chunks()):
                pulse_edges = chunk.meta['pulse_edges']
                spikes = detect_evoked_spikes(chunk, pulse_edges)
                spike_info.append({'pulse_n': chunk.meta['pulse_n'], 'pulse_start': pulse_edges[0], 'pulse_end': pulse_edges[1], 'spikes': spikes})
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


class PulseResponse(object):
    """Represents a chunk of postsynaptic recording taken during a presynaptic pulse stimulus.

    This class provides access to:
    - presynaptic stimulus waveform and metadata
    - presynaptic recording (if available)
    - postsynaptic recording
    - cell pair (if available)
    - postsynaptic cell
    """
    def __init__(self, recording=None, stim_pulse=None, pair=None, start_time=None, ex_qc_pass=None, in_qc_pass=None, post_tseries=None):
        self.recording = recording
        self.stim_pulse = stim_pulse
        self.pair = pair
        self.start_time = start_time
        self.ex_qc_pass = ex_qc_pass
        self.in_qc_pass = in_qc_pass
        self.post_tseries = post_tseries

    @property
    def pre_tseries(self):
        return self.stim_pulse.recorded_tseries

    @property
    def stim_tseries(self):
        return self.stim_pulse.stimulus_tseries
    

class StimPulse(object):
    """Represents a single stimiulus pulse intended to evoke a synaptic response.

    This class provides access to:
    - presynaptic stimulus waveform and metadata
    - presynaptic recording (if available)
    - presynaptic cell
    """
    def __init__(self, recording=None, pulse_number=None, onset_time=None, next_pulse_time=None, amplitude=None, duration=None, 
                 n_spikes=None, recorded_tseries=None):
        self.recording = recording
        self.pulse_number = pulse_number
        self.onset_time = onset_time
        self.next_pulse_time = next_pulse_time
        self.amplitude = amplitude
        self.duration = duration
        self.n_spikes = n_spikes
        self.recorded_tseries = recorded_tseries

    @property
    def recorded_tseries(self):
        if self._rec_tseries is None:
            self._rec_tseries = Trace(self.data, sample_rate=default_sample_rate, t0=self.data_start_time)
        return self._rec_tseries

    @property
    def stimulus_tseries(self):
        # can we work Stimulus objects into here, rather than generating manually??
        if self._stim_tseries is None:
            rec_ts = self.recorded_tseries  # todo: avoid loading full data for this step
            data = np.zeros(shape=rec_ts.shape)
            pstart = rec_ts.index_at(self.onset_time)
            pstop = rec_ts.index_at(self.onset_time + self.duration)
            data[pstart:pstop] = self.amplitude            
            self._stim_tseries = rec_ts.copy(data=data)
        return self._stim_tseries
