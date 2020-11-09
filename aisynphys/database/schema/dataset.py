from collections import OrderedDict
import numpy as np
from neuroanalysis.data import TSeries
from neuroanalysis.stimuli import Stimulus
from sqlalchemy.orm import relationship
from . import make_table
from .experiment import Experiment, Electrode, Pair, Cell
from . import default_sample_rate, sample_rate_str

__all__ = ['SyncRec', 'Recording', 'PatchClampRecording', 'MultiPatchProbe', 'TestPulse', 'StimPulse', 'StimSpike', 'PulseResponse', 'Baseline']


SyncRec = make_table(
    name='sync_rec',
    comment="A synchronous recording represents a \"sweep\" -- multiple recordings that were made simultaneously on different electrodes.",
    columns=[
        ('experiment_id', 'experiment.id', '', {'index': True}),
        ('ext_id', 'object', 'External ID of the SyncRecording'),
        ('temperature', 'float', 'Bath temperature during this recording'),
    ]
)

Experiment.sync_recs = relationship(SyncRec, order_by=SyncRec.id, back_populates="experiment", cascade='save-update,merge,delete', single_parent=True)
SyncRec.experiment = relationship(Experiment, back_populates='sync_recs')


class RecordingBase:
    def __repr__(self):
        sr_id = self.sync_rec.ext_id
        ex_id = self.sync_rec.experiment.ext_id
        return "<%s %s:%s.%s %s>" % (self.__class__.__name__, ex_id, sr_id, self.electrode.ext_id, self.stim_name)

    @property
    def stimulus(self):
        """An instance of neuroanalysis.stimuli.Stimulus describing the stimulus protocol used during
        this recording, or None if no stimulus information was recorded. 
        """
        stim = self.meta.get('stimulus', None)
        if stim is None:
            return None
        return Stimulus.load(stim)


Recording = make_table(
    name='recording',
    base=RecordingBase,
    comment= "A recording represents a single contiguous sweep recorded from a single electrode.",
    columns=[
        ('sync_rec_id', 'sync_rec.id', 'References the synchronous recording to which this recording belongs.', {'index': True}),
        ('electrode_id', 'electrode.id', 'Identifies the electrode that generated this recording', {'index': True}),
        ('start_time', 'datetime', 'The clock time at the start of this recording'),
        ('sample_rate', 'int', 'Sample rate for this recording'),
        ('device_name', 'str', 'Name of the device that generated this recording'),
        ('stim_name', 'str', 'The name of the stimulus protocol used in this recording, if any'),
        ('stim_meta', 'object', 'A data structure describing the stimulus protocol'),
    ]
)

SyncRec.recordings = relationship(Recording, order_by=Recording.id, back_populates="sync_rec", cascade='save-update,merge,delete', single_parent=True)
Recording.sync_rec = relationship(SyncRec, back_populates="recordings")

Electrode.recordings = relationship(Recording, back_populates="electrode", cascade='save-update,merge,delete', single_parent=True)
Recording.electrode = relationship(Electrode, back_populates="recordings")


PatchClampRecording = make_table(
    name='patch_clamp_recording',
    comment="Extra data for recordings made with a patch clamp amplifier",
    columns=[
        ('recording_id', 'recording.id', '', {'index': True, 'unique': True}),
        ('clamp_mode', 'str', 'The mode used by the patch clamp amplifier: "ic" or "vc"', {'index': True}),
        ('patch_mode', 'str', "The state of the membrane patch. E.g. 'whole cell', 'cell attached', 'loose seal', 'bath', 'inside out', 'outside out'"),
        ('baseline_potential', 'float', 'Median steady-state potential (recorded for IC or commanded for VC) during the recording'),
        ('baseline_current', 'float', 'Median steady-state current (recorded for VC or commanded for IC) during the recording'),
        ('baseline_rms_noise', 'float', 'RMS noise of the steady-state part of the recording'),
        ('nearest_test_pulse_id', 'test_pulse.id', 'ID of the test pulse that was recorded closest to this recording (and possibly embedded within the recording)'),
        ('qc_pass', 'bool', 'Indicates whether this recording passes a minimal ephys QC', {'index': True}),
    ]
)

Recording.patch_clamp_recording = relationship(PatchClampRecording, back_populates="recording", cascade='save-update,merge,delete', single_parent=True, uselist=False)
PatchClampRecording.recording = relationship(Recording, back_populates="patch_clamp_recording", single_parent=True)


MultiPatchProbe = make_table(
    name='multi_patch_probe',
    comment="Extra data for multipatch recordings intended to test synaptic dynamics.",
    columns=[
        ('patch_clamp_recording_id', 'patch_clamp_recording.id', '', {'index': True, 'unique': True}),
        ('induction_frequency', 'float', 'The induction frequency (Hz) of presynaptic pulses', {'index': True}),
        ('recovery_delay', 'float', 'The recovery delay (s) inserted between presynaptic pulses', {'index': True}),
        ('n_spikes_evoked', 'int', 'The number of presynaptic spikes evoked'),
    ]
)

PatchClampRecording.multi_patch_probe = relationship(MultiPatchProbe, back_populates="patch_clamp_recording", cascade='save-update,merge,delete', single_parent=True, uselist=False)
MultiPatchProbe.patch_clamp_recording = relationship(PatchClampRecording, back_populates="multi_patch_probe")


TestPulse = make_table(
    name='test_pulse',
    comment= "A short, usually hyperpolarizing pulse used to test the resistance of pipette, cell access, or cell membrane.",
    columns=[
        ('electrode_id', 'electrode.id', 'ID of the electrode on which this test pulse was recorded.', {'index': True}), 
        ('recording_id', 'recording.id', 'ID of the recording that contains this test pulse, if any.', {'index': True}),
        ('start_index', 'int'),
        ('stop_index', 'int'),
        ('baseline_current', 'float'),
        ('baseline_potential', 'float'),
        ('access_resistance', 'float'),
        ('input_resistance', 'float'),
        ('capacitance', 'float'),
        ('time_constant', 'float'),
    ]
)

Electrode.test_pulses = relationship(TestPulse, back_populates='electrode', cascade='save-update,merge,delete', single_parent=True)
TestPulse.electrode = relationship(Electrode, back_populates="test_pulses")
Recording.test_pulses = relationship(TestPulse, back_populates='recording', cascade='save-update,merge,delete', single_parent=True)
TestPulse.recording = relationship(Recording, back_populates="test_pulses")
PatchClampRecording.nearest_test_pulse = relationship(TestPulse, single_parent=True, foreign_keys=[PatchClampRecording.nearest_test_pulse_id])


class StimPulseBase(object):
    def _init_on_load(self):
        self._rec_tseries = None
        self._stim_tseries = None
    
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

   
StimPulse = make_table(
    name='stim_pulse',
    base=StimPulseBase,
    comment= "A pulse stimulus intended to evoke an action potential",
    columns=[
        ('recording_id', 'recording.id', '', {'index': True}),
        ('pulse_number', 'int', 'The ordinal position of this pulse within a train of pulses.', {'index': True}),
        ('cell_id', 'cell.id', 'Cell that was targeted by this stimulus, if any.', {'index':True}),
        ('onset_time', 'float', 'The starting time of the pulse, relative to the beginning of the recording'),
        ('amplitude', 'float', 'Amplitude of the presynaptic pulse'),
        ('duration', 'float', 'Length of the pulse in seconds'),
        ('n_spikes', 'int', 'Number of spikes evoked by this pulse'),
        ('first_spike_time', 'float', 'Time of the first spike evoked by this pulse, measured from the beginning of the recording until the max slope of the spike rising phase.'),
        # ('first_spike', 'stim_spike.id', 'The ID of the first spike evoked by this pulse'),
        ('data', 'array', 'Numpy array of presynaptic recording sampled at '+sample_rate_str, {'deferred': True}),
        ('data_start_time', 'float', "Starting time of the data chunk, relative to the beginning of the recording"),       
        ('position', 'object', '3D location of this stimulation in the arbitrary coordinate system of the experiment'),
        ('qc_pass', 'bool', 'Indicates whether this stimulation passed qc.'),
        ('previous_pulse_dt', 'float', 'Time elapsed since the last stimulus in the same cell', {'index': True}),
    ]
)

Recording.stim_pulses = relationship(StimPulse, back_populates="recording", cascade='save-update,merge,delete', single_parent=True)
StimPulse.recording = relationship(Recording, back_populates="stim_pulses")
Cell.stim_pulses = relationship(StimPulse, back_populates="cell", cascade='save-update,merge,delete', single_parent=True)
StimPulse.cell = relationship(Cell, back_populates='stim_pulses')


StimSpike = make_table(
    name='stim_spike',
    comment= "An action potential evoked by a stimulus pulse. Note that some metrics may be omitted if they could not be determined accurately.",
    columns=[
        ('stim_pulse_id', 'stim_pulse.id', '', {'index': True}),
        ('onset_time', 'float', "The time of the earliest detectable effect of the spike."),
        ('max_slope_time', 'float', "The time of the max slope of the spike, relative to the beginning of the recording."),
        ('max_slope', 'float', 'Maximum slope of the presynaptic spike'),
        ('peak_time', 'float', "The time of the peak of the spike, relative to the beginning of the recording."),
        ('peak_diff', 'float', 'Amplitude of the spike peak, relative to baseline'),
        ('peak_value', 'float', 'Absolute value of the spike peak'),
    ]
)

StimSpike.stim_pulse = relationship(StimPulse, back_populates="spikes")
StimPulse.spikes = relationship(StimSpike, back_populates="stim_pulse", single_parent=True)


class PulseResponseBase(object):
    def _init_on_load(self):
        self._post_tseries = None
    
    @property
    def post_tseries(self):
        if self._post_tseries is None:
            self._post_tseries = TSeries(self.data, sample_rate=default_sample_rate, t0=self.data_start_time)
        return self._post_tseries

    @property
    def pre_tseries(self):
        return self.stim_pulse.recorded_tseries

    @property
    def stim_tseries(self):
        return self.stim_pulse.stimulus_tseries
        
    @property
    def baseline_tseries(self):
        bl = self.baseline
        if bl is None:
            return None
        return TSeries(bl.data, sample_rate=default_sample_rate, t0=bl.data_start_time)

    def get_tseries(self, ts_type, align_to):
        """Return the pre-, post-, or stimulus TSeries, time aligned to either the spike or the stimulus onset.
            
        If spike alignment is requested but no spike time is available, then return None.

        If any alignment is requested, this method also adds keys to the returned ``tseries.meta``: 
        'pulse_start', 'pulse_stop', and 'spike_time', giving the times of these events relative to the aligned timebase.
        
        Parameters
        ----------
        ts_type : str
            One of "pre", "post", "stim", or "baseline"
        align_to : str | None
            One of "spike" or "pulse"        
        """
        assert ts_type in ('pre', 'post', 'stim', 'baseline'), "ts_type must be 'pre', 'post', 'stim', or 'baseline'"
        assert align_to in ('spike', 'pulse', None), "align_to must be 'spike', 'stim', or None"
        
        pulse_time = self.stim_pulse.onset_time
        spike_time = self.stim_pulse.first_spike_time
        
        ts = getattr(self, ts_type+"_tseries")
        if ts is None or align_to == None:
            return ts

        if ts_type == 'baseline':
            # if pulse/spike aligning, we set up timing to be exactly
            # the same as the post recording
            ts = ts.copy(t0=self.data_start_time)
        
        if align_to == 'spike':
            align_time = self.stim_pulse.first_spike_time
            if align_time is None:
                return None
        elif align_to == 'pulse':
            align_time = self.stim_pulse.onset_time
            
        ts = ts.copy(t0=ts.t0 - align_time)
        ts.meta.update({
            'pulse_start': pulse_time - align_time,
            'pulse_stop': pulse_time - align_time + self.stim_pulse.duration,
            'spike_time': None if spike_time is None else spike_time - align_time,
        })
        return ts
        

PulseResponse = make_table(
    name='pulse_response',
    base=PulseResponseBase,
    comment="A chunk of postsynaptic recording taken during a presynaptic pulse stimulus",
    columns=[
        ('recording_id', 'recording.id', 'The full recording from which this pulse was extracted', {'index': True}),
        ('stim_pulse_id', 'stim_pulse.id', 'The presynaptic pulse', {'index': True}),
        ('pair_id', 'pair.id', 'The pre-post cell pair involved in this pulse response', {'index': True}),
        ('baseline_id', 'baseline.id', 'A random baseline snippet matched from the same recording.', {'index': True}),
        ('data', 'array', 'numpy array of response data sampled at '+sample_rate_str, {'deferred': True}),
        ('data_start_time', 'float', 'Starting time of this chunk of the recording in seconds, relative to the beginning of the recording'),
        ('ex_qc_pass', 'bool', 'Indicates whether this recording snippet passes QC for excitatory synapse probing', {'index': True}),
        ('in_qc_pass', 'bool', 'Indicates whether this recording snippet passes QC for inhibitory synapse probing', {'index': True}),
    ]
)

StimPulse.pulse_response = relationship(PulseResponse, back_populates="stim_pulse", cascade='save-update,merge,delete', single_parent=True)
PulseResponse.stim_pulse = relationship(StimPulse)
PulseResponse.recording = relationship(Recording)
Pair.pulse_responses = relationship(PulseResponse, back_populates='pair', single_parent=True, order_by=PulseResponse.id)
PulseResponse.pair = relationship(Pair, back_populates='pulse_responses')


Baseline = make_table(
    name='baseline',
    comment="A snippet of baseline data used for comparison to pulse_response records",
    columns=[
        ('recording_id', 'recording.id', 'The recording from which this baseline snippet was extracted.', {'index': True}),
        ('data', 'array', 'numpy array of baseline data sampled at '+sample_rate_str, {'deferred': True}),
        ('data_start_time', 'float', "Starting time of this chunk of the recording in seconds, relative to the beginning of the recording"),
        ('mode', 'float', 'most common value in the baseline snippet'),
        ('ex_qc_pass', 'bool', 'Indicates whether this recording snippet passes QC for excitatory synapse probing'),
        ('in_qc_pass', 'bool', 'Indicates whether this recording snippet passes QC for inhibitory synapse probing'),
    ]
)

Recording.baselines = relationship(Baseline, back_populates="recording", cascade='save-update,merge,delete', single_parent=True)
Baseline.recording = relationship(Recording, back_populates="baselines")
PulseResponse.baseline = relationship(Baseline, back_populates="pulse_responses", cascade='save-update,merge,delete', uselist=False)
Baseline.pulse_responses = relationship(PulseResponse, back_populates="baseline", uselist=True)
