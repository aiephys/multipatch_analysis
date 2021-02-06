import sys
import numpy as np

from neuroanalysis.miesnwb import MiesNwb, MiesSyncRecording, MiesRecording
from neuroanalysis.stimuli import find_square_pulses
from neuroanalysis.spike_detection import detect_evoked_spikes
from neuroanalysis.data import TSeries, TSeriesList
from neuroanalysis.baseline import float_mode
from neuroanalysis.analyzers.analyzer import Analyzer
from neuroanalysis.analyzers.stim_pulse import PatchClampStimPulseAnalyzer
from neuroanalysis.analyzers.baseline import BaselineDistributor

from .. import qc


class MultiPatchDataset(MiesNwb):
    """Extension of neuroanalysis data abstraction layer to include
    multipatch-specific metadata.
    """
    def create_sync_recording(self, sweep_id):
        return MultiPatchSyncRecording(self, sweep_id)

        
class MultiPatchSyncRecording(MiesSyncRecording):
    def __init__(self, nwb, sweep_id):
        MiesSyncRecording.__init__(self, nwb, sweep_id)
        self._baseline_mask = None
        self._baseline_regions = None
        try:
            self.meta['temperature'] = self.recordings[0].meta['notebook']['Async AD 1: Bath Temperature']
        except Exception:
            self.meta['temperature'] = None
    
    def create_recording(self, sweep_id, ch):
        miesrec = MiesRecording(self, sweep_id, ch)
        stim = miesrec.meta['notebook']['Stim Wave Name'].lower()
        if any(substr in stim for substr in ['pulsetrain', 'recovery', 'pulstrn']):
            return MultiPatchProbe(miesrec)
        elif any(substr in stim for substr in ['mixedfs']):
            return MultiPatchMixedFreqTrain(miesrec)
        else:
            return MultiPatchRecording(miesrec)

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
                pa = PatchClampStimPulseAnalyzer.get(rec)
                pulses = pa.pulses()
                for pulse in pulses:
                    start = pri.index_at(pulse[0])
                    stop = pri.index_at(pulse[1])
                    mask[start:stop + settle_size] = True
            # mask is False in quiescent regions
            self._baseline_mask = mask

        if self._baseline_regions is None:
            mask = self._baseline_mask
            starts = list(np.argwhere(~mask[1:] & mask[:-1])[:,0])
            if not mask[0]:
                starts.insert(0, 0)
            stops = list(np.argwhere(mask[1:] & ~mask[:-1])[:,0])
            if not mask[-1]:
                stops.append(len(mask))
            
            baseline_inds = [r for r in zip(starts, stops) if r[1] > r[0]]
            self._baseline_regions = [(pri.time_at(i0), pri.time_at(i1)) for i0, i1 in baseline_inds]

        return self._baseline_regions


class MultiPatchRecording(MiesRecording):
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
        # ask the parent sweep for baseline regions in which no channels are active
        if self._base_regions is None:
            self._base_regions = self._parent_rec.parent.baseline_regions()
        return self._base_regions


class MultiPatchProbe(MultiPatchRecording):
    """A 12-pulse stimulus/response used to probe for synaptic connections.
    """
    def stim_params(self):
        """Return induction frequency and recovery delay.
        """
        psa = PatchClampStimPulseAnalyzer.get(self)

        pulses = [p[0] for p in psa.pulses(channel='command') if p[2] > 0]
        if len(pulses) < 2:
            return None, None
        ind_freq = np.round(1.0 / (pulses[1] - pulses[0]))
        rec_delay = np.round(np.diff(pulses).max(), 3)

        return ind_freq, rec_delay


class MultiPatchMixedFreqTrain(MultiPatchRecording):
    """A spike train composed from random interspike intevals.
    """
    pass


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
        
        Returns
        -------
        responses : list
            Each item in the list is a dict corresponding to a single presynaptic stimulus, with the following keys:
            * pulse_n: number of this stimulus within the sweep
            * pulse_start: start time of the stimulus pulse
            * pulse_end: end time of the stimulus pulse 
            * spikes: list of presynaptic spikes detected 
            * response: time slice of the postsynaptic recording
            * pre_rec: time slice of the presynaptic recording
        
        """
        # detect presynaptic spikes
        pulse_stim = PatchClampStimPulseAnalyzer.get(pre_rec)
        spikes = pulse_stim.evoked_spikes()
        # spikes looks like:
        #   [{'pulse_n', 'pulse_start', 'pulse_end', 'spikes': [...]}, ...]

        # Select ranges to extract from postsynaptic recording
        result = []
        for i,pulse in enumerate(spikes):
            pulse = pulse.copy()
            if len(pulse['spikes']) == 0:
                if require_spike:
                    continue
                spike = None
            else:
                spike = pulse['spikes'][0]
            
            if align_to == 'spike':
                # start recording window at the rising phase of the presynaptic spike
                pulse['rec_start'] = spike['max_slope_time'] - pre_pad
            elif align_to == 'pulse':
                # align to pulse onset
                pulse['rec_start'] = pulse['pulse_start'] - pre_pad
            pulse['rec_start'] = max(pre_rec['primary'].t0, pulse['rec_start'])
            
            # get times of nearby pulses
            prev_pulse = None if i == 0 else spikes[i-1]['pulse_end']
            this_pulse = pulse['pulse_start']
            next_pulse = None if i+1 >= len(spikes) else spikes[i+1]['pulse_start']

            max_stop = pulse['rec_start'] + 50e-3
            if next_pulse is not None:
                # truncate window early if there is another pulse
                pulse['rec_stop'] = min(max_stop, next_pulse)
            else:
                # otherwise, stop 50 ms later
                pulse['rec_stop'] = max_stop
            
            # Extract data from postsynaptic recording
            pulse['response'] = post_rec.time_slice(pulse['rec_start'], pulse['rec_stop'])
            assert len(pulse['response']['primary']) > 0

            # Extract presynaptic spike and stimulus command
            pulse['pre_rec'] = pre_rec.time_slice(pulse['rec_start'], pulse['rec_stop'])

            # Add minimal QC metrics for excitatory and inhibitory measurements
            pulse_window = [pulse['rec_start'], pulse['rec_stop']]
            n_spikes = 0 if spike is None else 1
            adj_pulse_times = []
            if prev_pulse is not None:
                adj_pulse_times.append(prev_pulse - this_pulse)
            if next_pulse is not None:
                adj_pulse_times.append(next_pulse - this_pulse)
            pulse['ex_qc_pass'], pulse['in_qc_pass'], pulse['qc_failures'] = qc.pulse_response_qc_pass(post_rec=post_rec, window=pulse_window, n_spikes=n_spikes, adjacent_pulses=adj_pulse_times)

            result.append(pulse)
        
        return result

    def get_pulse_response(self, pre_rec, post_rec, first_pulse=0, last_pulse=-1):
        pulse_stim = PatchClampStimPulseAnalyzer.get(pre_rec)
        spikes = pulse_stim.evoked_spikes()
        
        if len(spikes) < 10:
            return None
        
        start = spikes[first_pulse]['pulse_start'] - 20e-3
        stop = spikes[last_pulse]['pulse_start'] + 50e-3
        return post_rec['primary'].time_slice(start, stop)

    def get_train_response(self, pre_rec, post_rec, start_pulse, stop_pulse, padding=(-10e-3, 50e-3)):
        """Return the part of the post-synaptic recording during a range of pulses,
        along with a baseline chunk
        """
        pulse_stim = PatchClampStimPulseAnalyzer.get(pre_rec)
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

    def find_artifacts(self, rec, freq=None, pos_threshold=-10e-3, neg_threshold=-100e-3):
        """Find contaminating artifacts in recording such as spontaneous spikes or electrical noise, return True
        state if anything is found"""

        artifacts = False
        # find anything that crosses a max threshold of -10mV
        if np.any(rec >= pos_threshold):
            artifacts = True
        if np.any(rec <= neg_threshold):
            artifacts = True
        return artifacts


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
    

class PulseResponseList(object):
    """A list of pulse responses with methods for time-aligning and baseline
    subtracting recordings.

    Parameters
    ----------
    prs : list
        A list of PulseResponse instances
    """
    def __init__(self, prs):
        self.prs = prs

    def __len__(self):
        return len(self.prs)

    def __getitem__(self, item):
        return self.prs[item]

    def __iter__(self):
        for pr in self.prs:
            yield pr

    def post_tseries(self, align=None, bsub=False, bsub_win=5e-3, alignment_failure_mode='ignore'):
        """Return a TSeriesList of all postsynaptic recordings.
        """
        return self._get_tserieslist('post_tseries', align, bsub, bsub_win, alignment_failure_mode)

    def pre_tseries(self, align=None, bsub=False, bsub_win=5e-3, alignment_failure_mode='ignore'):
        """Return a TSeriesList of all presynaptic recordings.
        """
        return self._get_tserieslist('pre_tseries', align, bsub, bsub_win, alignment_failure_mode)

    def _get_tserieslist(self, ts_name, align, bsub, bsub_win=5e-3, alignment_failure_mode='ignore'):
        tsl = []
        if align is not None and alignment_failure_mode == 'average':
            if align == 'spike':
                average_align_t = np.mean([p.stim_pulse.first_spike_time for p in self.prs if p.stim_pulse.first_spike_time is not None])
            elif align == 'peak':
                average_align_t = np.mean([p.stim_pulse.spikes[0].peak_time for p in self.prs if p.stim_pulse.n_spikes==1 and p.stim_pulse.spikes[0].peak_time is not None])
            elif align == 'stim':
                average_align_t = np.mean([p.stim_pulse.onset_time for p in self.prs if p.stim_pulse.onset_time is not None])
            else:
                raise ValueError("align must be None, 'spike', 'peak', or 'pulse'.")
        
        for pr in self.prs:   
            ts = getattr(pr, ts_name)
            stim_time = pr.stim_pulse.onset_time

            if bsub is True:
                start_time = max(ts.t0, stim_time-bsub_win)
                baseline_data = ts.time_slice(start_time, stim_time).data
                if len(baseline_data) == 0:
                    baseline = ts.data[0]
                else:
                    baseline = float_mode(baseline_data)
                ts = ts - baseline

            if align is not None:
                if align == 'spike':
                # first_spike_time is the max dv/dt of the spike
                    align_t = pr.stim_pulse.first_spike_time
                elif align == 'pulse':
                    align_t = stim_time
                elif align == 'peak':
                    # peak of the first spike
                    align_t = pr.stim_pulse.spikes[0].peak_time if pr.stim_pulse.n_spikes==1 else None
                else:
                    raise ValueError("align must be None, 'spike', 'peak', or 'pulse'.")
                
                if align_t is None:
                    if alignment_failure_mode == 'ignore':
                        # ignore PRs with no known timing
                        continue
                    elif alignment_failure_mode == 'average':
                        align_t = average_align_t
                        if np.isnan(align_t):
                            raise Exception("average %s time is None, try another mode" % align)
                    elif alignment_failure_mode == 'raise':
                        raise Exception("%s time is not available for pulse %s and can't be aligned" % (align, pr))
                
                ts = ts.copy(t0=ts.t0 - align_t)

            tsl.append(ts)
        return TSeriesList(tsl)


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
            self._rec_tseries = TSeries(self.data, sample_rate=default_sample_rate, t0=self.data_start_time)
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
