from .database import TableGroup, _generate_mapping, _sample_rate_str
from sqlalchemy.orm import relationship, deferred, sessionmaker, aliased
from .experiment import Experiment, Electrode, Pair


class DatasetTableGroup(TableGroup):
    """Defines tables that hold data imported from a dataset, such as recordings, stimuli, pulse responses, etc.
    """
    schemas = {
        'sync_rec': [
            """A synchronous recording represents a "sweep" -- multiple recordings that were made simultaneously
            on different electrodes.""",
            ('experiment_id', 'experiment.id', '', {'index': True}),
            ('ext_id', 'object', 'External ID of the SyncRecording'),
            ('temperature', 'float', 'Bath temperature during this recording'),
        ],
        'recording': [
            """A recording represents a single contiguous sweep recorded from a single electrode. 
            """,
            ('sync_rec_id', 'sync_rec.id', 'References the synchronous recording to which this recording belongs.', {'index': True}),
            ('electrode_id', 'electrode.id', 'Identifies the electrode that generated this recording', {'index': True}),
            ('start_time', 'datetime', 'The clock time at the start of this recording'),
            ('sample_rate', 'int', 'Sample rate for this recording'),
        ],
        'patch_clamp_recording': [
            "Extra data for recordings made with a patch clamp amplifier",
            ('recording_id', 'recording.id', '', {'index': True, 'unique': True}),
            ('clamp_mode', 'str', 'The mode used by the patch clamp amplifier: "ic" or "vc"', {'index': True}),
            ('patch_mode', 'str', "The state of the membrane patch. E.g. 'whole cell', 'cell attached', 'loose seal', 'bath', 'inside out', 'outside out'"),
            ('stim_name', 'str', "The name of the stimulus protocol"),
            ('baseline_potential', 'float', 'Median steady-state potential (recorded for IC or commanded for VC) during the recording'),
            ('baseline_current', 'float', 'Median steady-state current (recorded for VC or commanded for IC) during the recording'),
            ('baseline_rms_noise', 'float', 'RMS noise of the steady-state part of the recording'),
            ('nearest_test_pulse_id', 'test_pulse.id', 'ID of the test pulse that was recorded closest to this recording (and possibly embedded within the recording)'),
            ('qc_pass', 'bool', 'Indicates whether this recording passes a minimal ephys QC'),
        ],
        'multi_patch_probe': [
            "Extra data for multipatch recordings intended to test synaptic dynamics.",
            ('patch_clamp_recording_id', 'patch_clamp_recording.id', '', {'index': True, 'unique': True}),
            ('induction_frequency', 'float', 'The induction frequency (Hz) of presynaptic pulses', {'index': True}),
            ('recovery_delay', 'float', 'The recovery delay (s) inserted between presynaptic pulses', {'index': True}),
            ('n_spikes_evoked', 'int', 'The number of presynaptic spikes evoked'),
        ],
        'test_pulse': [
            """A short, usually hyperpolarizing pulse used to test the resistance of pipette, cell access, or cell membrane.
            """,
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
        ],
        'stim_pulse': [
            "A pulse stimulus intended to evoke an action potential",
            ('recording_id', 'recording.id', '', {'index': True}),
            ('pulse_number', 'int', 'The ordinal position of this pulse within a train of pulses.', {'index': True}),
            ('onset_time', 'float', 'The starting time of the pulse, relative to the beginning of the recording'),
            ('next_pulse_time', 'float', 'Time of the next pulse on any channel in the sync rec'),
            ('amplitude', 'float', 'Amplitude of the presynaptic pulse'),
            ('duration', 'float', 'Length of the pulse in seconds'),
            ('n_spikes', 'int', 'Number of spikes evoked by this pulse'),
            # ('first_spike', 'stim_spike.id', 'The ID of the first spike evoked by this pulse'),
            ('data', 'array', 'Numpy array of presynaptic recording sampled at '+_sample_rate_str, {'deferred': True}),
            ('data_start_time', 'float', "Starting time of the data chunk, relative to the beginning of the recording"),
        ],
        'stim_spike': [
            "An action potential evoked by a stimulus pulse",
            ('stim_pulse_id', 'stim_pulse.id', '', {'index': True}),
            ('peak_time', 'float', "The time of the peak of the spike, relative to the beginning of the recording."),
            ('peak_diff', 'float', 'Amplitude of the spike peak, relative to baseline'),
            ('peak_val', 'float', 'Absolute value of the spike peak'),
            ('max_dvdt_time', 'float', "The time of the max dv/dt of the spike, relative to the beginning of the recording."),
            ('max_dvdt', 'float', 'Maximum slope of the presynaptic spike'),
        ],
        'baseline': [
            "A snippet of baseline data, matched to a postsynaptic recording",
            ('recording_id', 'recording.id', 'The recording from which this baseline snippet was extracted.', {'index': True}),
            ('start_time', 'float', "Starting time of this chunk of the recording in seconds, relative to the beginning of the recording"),
            ('data', 'array', 'numpy array of baseline data sampled at '+_sample_rate_str, {'deferred': True}),
            ('mode', 'float', 'most common value in the baseline snippet'),
            ('ex_qc_pass', 'bool', 'Indicates whether this recording snippet passes QC for excitatory synapse probing'),
            ('in_qc_pass', 'bool', 'Indicates whether this recording snippet passes QC for inhibitory synapse probing'),
        ],
        'pulse_response': [
            "A chunk of postsynaptic recording taken during a presynaptic pulse stimulus",
            ('recording_id', 'recording.id', 'The full recording from which this pulse was extracted', {'index': True}),
            ('stim_pulse_id', 'stim_pulse.id', 'The presynaptic pulse', {'index': True}),
            ('pair_id', 'pair.id', 'The pre-post cell pair involved in this pulse response', {'index': True}),
            ('start_time', 'float', 'Starting time of this chunk of the recording in seconds, relative to the beginning of the recording'),
            ('data', 'array', 'numpy array of response data sampled at '+_sample_rate_str, {'deferred': True}),
            ('ex_qc_pass', 'bool', 'Indicates whether this recording snippet passes QC for excitatory synapse probing'),
            ('in_qc_pass', 'bool', 'Indicates whether this recording snippet passes QC for inhibitory synapse probing'),
        ],
    }

    def create_mappings(self):
        TableGroup.create_mappings(self)
        SyncRec = self['sync_rec']
        Recording = self['recording']
        PatchClampRecording = self['patch_clamp_recording']
        MultiPatchProbe = self['multi_patch_probe']
        TestPulse = self['test_pulse']
        StimPulse = self['stim_pulse']
        StimSpike = self['stim_spike']
        PulseResponse = self['pulse_response']
        Baseline = self['baseline']

        # Set up relationships
        Experiment.sync_recs = relationship(SyncRec, order_by=SyncRec.id, back_populates="experiment", cascade='save-update,merge,delete', single_parent=True)
        SyncRec.experiment = relationship(Experiment, back_populates='sync_recs')

        SyncRec.recordings = relationship(Recording, order_by=Recording.id, back_populates="sync_rec", cascade='save-update,merge,delete', single_parent=True)
        Recording.sync_rec = relationship(SyncRec, back_populates="recordings")

        Electrode.recordings = relationship(Recording, back_populates="electrode", cascade='save-update,merge,delete', single_parent=True)
        Recording.electrode = relationship(Electrode, back_populates="recordings")

        Recording.patch_clamp_recording = relationship(PatchClampRecording, back_populates="recording", cascade='save-update,merge,delete', single_parent=True, uselist=False)
        PatchClampRecording.recording = relationship(Recording, back_populates="patch_clamp_recording", single_parent=True)

        PatchClampRecording.multi_patch_probe = relationship(MultiPatchProbe, back_populates="patch_clamp_recording", cascade='save-update,merge,delete', single_parent=True)
        MultiPatchProbe.patch_clamp_recording = relationship(PatchClampRecording, back_populates="multi_patch_probe")

        Electrode.test_pulses = relationship(TestPulse, back_populates='electrode', cascade='save-update,merge,delete', single_parent=True)
        TestPulse.electrode = relationship(Electrode, back_populates="test_pulses")
        Recording.test_pulses = relationship(TestPulse, back_populates='recording', single_parent=True)
        TestPulse.recording = relationship(Recording, back_populates="test_pulses")

        PatchClampRecording.nearest_test_pulse = relationship(TestPulse, single_parent=True, foreign_keys=[PatchClampRecording.nearest_test_pulse_id])
        #TestPulse.patch_clamp_recording = relationship(PatchClampRecording)

        Recording.stim_pulses = relationship(StimPulse, back_populates="recording", cascade='save-update,merge,delete', single_parent=True)
        StimPulse.recording = relationship(Recording, back_populates="stim_pulses")

        StimSpike.pulse = relationship(StimPulse, back_populates="spikes")
        StimPulse.spikes = relationship(StimSpike, back_populates="pulse", single_parent=True)
        StimPulse.pulse_response = relationship(PulseResponse, back_populates="stim_pulse", cascade='save-update,merge,delete', single_parent=True)
        
        Recording.baselines = relationship(Baseline, back_populates="recording", cascade='save-update,merge,delete', single_parent=True)
        Baseline.recording = relationship(Recording, back_populates="baselines")

        PulseResponse.recording = relationship(Recording)
        PulseResponse.stim_pulse = relationship(StimPulse)
        Pair.pulse_responses = relationship(PulseResponse, back_populates='pair', single_parent=True)
        PulseResponse.pair = relationship(Pair, back_populates='pulse_responses')


dataset_tables = DatasetTableGroup()


def init_tables():
    global SyncRec, Recording, PatchClampRecording, MultiPatchProbe, TestPulse, StimPulse, StimSpike, PulseResponse, Baseline
    dataset_tables.create_tables()

    SyncRec = dataset_tables['sync_rec']
    Recording = dataset_tables['recording']
    PatchClampRecording = dataset_tables['patch_clamp_recording']
    MultiPatchProbe = dataset_tables['multi_patch_probe']
    TestPulse = dataset_tables['test_pulse']
    StimPulse = dataset_tables['stim_pulse']
    StimSpike = dataset_tables['stim_spike']
    PulseResponse = dataset_tables['pulse_response']
    Baseline = dataset_tables['baseline']


# create tables in database and add global variables for ORM classes
init_tables()
