"""
Accumulate all experiment data into a set of linked tables.
"""
import io
import numpy as np

import sqlalchemy
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Boolean, Float, Date, DateTime, LargeBinary, ForeignKey
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship
from sqlalchemy.types import TypeDecorator
from sqlalchemy.orm import sessionmaker, aliased
from sqlalchemy.sql.expression import func

from .. import config

default_sample_rate = 20000


_sample_rate_str = '%dkHz' % (default_sample_rate // 1000)

table_schemas = {
    'slice': [
        "All brain slices on which an experiment was attempted.",
        ('acq_timestamp', 'datetime', 'Creation timestamp for slice data acquisition folder.', {'unique': True}),
        ('species', 'str', 'Human | mouse (from LIMS)'),
        ('age', 'int', 'Specimen age (in days) at time of dissection (from LIMS)'),
        ('sex', 'str', 'Specimen sex ("M", "F", or "unknown"; from LIMS)'),
        ('weight', 'str', 'Specimen weight (from LIMS)'),
        ('genotype', 'str', 'Specimen donor genotype (from LIMS)'),
        ('orientation', 'str', 'Orientation of the slice plane (eg "sagittal"; from LIMS specimen name)'),
        ('surface', 'str', 'The surface of the slice exposed during the experiment (eg "left"; from LIMS specimen name)'),
        ('hemisphere', 'str', 'The brain hemisphere from which the slice originated. (from LIMS specimen name)'),
        ('quality', 'int', 'Experimenter subjective slice quality assessment (0-5)'),
        ('slice_time', 'datetime', 'Time when this specimen was sliced'),
        ('slice_conditions', 'object', 'JSON containing solutions, perfusion, incubation time, etc.'),
        ('lims_specimen_name', 'str', 'Name of LIMS "slice" specimen'),
        ('storage_path', 'str', 'Location of data within server or cache storage'),
        ('submission_data', 'object'),          # structure generated for original submission
    ],
    'experiment': [
        "A group of cells patched simultaneously in the same slice.",
        ('original_path', 'str', 'Original location of raw data on rig.'),
        ('storage_path', 'str', 'Location of data within server or cache storage.'),
        ('ephys_file', 'str', 'Name of ephys NWB file relative to storage_path.'),
        ('rig_name', 'str', 'Identifier for the rig that generated these results.'),
        ('acq_timestamp', 'datetime', 'Creation timestamp for site data acquisition folder.', {'unique': True, 'index': True}),
        ('slice_id', 'slice.id', 'ID of the slice used for this experiment'),
        ('target_region', 'str', 'The intended brain region for this experiment'),
        ('internal', 'str', 'The name of the internal solution used in this experiment. '
                            'The solution should be described in the pycsf database.'),
        ('acsf', 'str', 'The name of the ACSF solution used in this experiment. '
                        'The solution should be described in the pycsf database.'),
        ('target_temperature', 'float', 'The intended temperature of the experiment (but actual recording temperature is stored elsewhere)'),
        ('date', 'datetime', 'The date of this experiment'),
        ('lims_specimen_id', 'int', 'ID of LIMS "CellCluster" specimen.'),
        ('submission_data', 'object', 'structure generated for original submission.'),
        ('lims_trigger_id', 'int', 'ID used to query status of LIMS upload.'),
        ('connectivity_analysis_complete', 'bool'),
        ('kinetics_analysis_complete', 'bool'),
    ],
    'electrode': [
        "Each electrode records a patch attempt, whether or not it resulted in a "
        "successful cell recording.",
        ('expt_id', 'experiment.id', '', {'index': True}),
        ('patch_status', 'str', 'no seal, low seal, GOhm seal, tech fail, ...'),
        ('start_time', 'datetime'),
        ('stop_time', 'datetime'),
        ('device_id', 'int', 'External identifier for the device attached to this electrode (usually the MIES A/D channel)'),
        ('initial_resistance', 'float'),
        ('initial_current', 'float'),
        ('pipette_offset', 'float'),
        ('final_resistance', 'float'),
        ('final_current', 'float'),
        ('notes', 'str'),
        ('ext_id', 'int', 'Electrode ID (usually 1-8) referenced in external metadata records'),
    ],
    'cell': [
        "Each row represents a single patched cell.",
        ('electrode_id', 'electrode.id'),
        ('cre_type', 'str', 'Comma-separated list of cre drivers apparently expressed by this cell'),
        ('target_layer', 'str', 'The intended cortical layer for this cell (used as a placeholder until the actual layer call is made)'),
        ('is_excitatory', 'bool', 'True if the cell is determined to be excitatory by synaptic current, cre type, or morphology'),
        ('synapse_sign', 'int', 'The sign of synaptic currents produced by this cell: excitatory=+1, inhibitory=-1, mixed=0'),
        ('patch_start', 'float'),
        ('patch_stop', 'float'),
        ('seal_resistance', 'float', 'The seal resistance recorded for this cell immediately before membrane rupture'),
        ('has_biocytin', 'bool', 'If true, then the soma was seen to be darkly stained with biocytin (this indicates a good reseal, but does may not indicate a high-quality fill)'),
        ('has_dye_fill', 'bool', 'Indicates whether the cell was filled with fluorescent dye during the experiment'),
        ('pass_qc', 'bool'),
        ('pass_spike_qc', 'bool'),
        ('depth', 'float', 'Depth of the cell (in m) from the cut surface of the slice.'),
        ('position', 'object', '3D location of this cell in the arbitrary coordinate system of the experiment'),
        ('ext_id', 'int', 'Cell ID (usually 1-8) referenced in external metadata records'),
    ],
    'pair': [
        "All possible putative synaptic connections",
        ('expt_id', 'experiment.id', '', {'index': True}),
        ('pre_cell_id', 'cell.id', 'ID of the presynaptic cell', {'index': True}),
        ('post_cell_id', 'cell.id', 'ID of the postsynaptic cell', {'index': True}),
        ('synapse', 'bool', 'Whether the experimenter thinks there is a synapse', {'index': True}),
        ('electrical', 'bool', 'Whether the experimenter thinks there is a gap junction', {'index': True}),
        ('crosstalk_artifact', 'float', 'Amplitude of crosstalk artifact measured in current clamp'),
        ('n_ex_test_spikes', 'int', 'Number of QC-passed spike-responses recorded for this pair at excitatory holding potential'),
        ('n_in_test_spikes', 'int', 'Number of QC-passed spike-responses recorded for this pair at inhibitory holding potential'),
        ('synapse_sign', 'int', 'Sign of synaptic current amplitude (+1 for excitatory, -1 for inhibitory'),
        ('distance', 'float', 'Distance between somas (in m)'),
    ],
    'sync_rec': [
        """A synchronous recording represents a "sweep" -- multiple recordings that were made simultaneously
        on different electrodes.""",
        ('experiment_id', 'experiment.id', '', {'index': True}),
        ('sync_rec_key', 'object', 'External ID of the SyncRecording'),
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
        ('recording_id', 'recording.id', '', {'index': True}),
        ('electrode_id', 'electrode.id', 'References the patch electrode that was used during this recording', {'index': True}),
        ('clamp_mode', 'str', 'The mode used by the patch clamp amplifier: "ic" or "vc"', {'index': True}),
        ('patch_mode', 'str', "The state of the membrane patch. E.g. 'whole cell', 'cell attached', 'loose seal', 'bath', 'inside out', 'outside out'"),
        ('stim_name', 'object', "The name of the stimulus protocol"),
        ('baseline_potential', 'float', 'Median steady-state potential (recorded for IC or commanded for VC) during the recording'),
        ('baseline_current', 'float', 'Median steady-state current (recorded for VC or commanded for IC) during the recording'),
        ('baseline_rms_noise', 'float', 'RMS noise of the steady-state part of the recording'),
        ('nearest_test_pulse_id', 'test_pulse.id', 'ID of the test pulse that was recorded closest to this recording (and possibly embedded within the recording)'),
    ],
    'multi_patch_probe': [
        "Extra data for multipatch recordings intended to test synaptic dynamics.",
        ('patch_clamp_recording_id', 'patch_clamp_recording.id', '', {'index': True}),
        ('induction_frequency', 'float', 'The induction frequency (Hz) of presynaptic pulses', {'index': True}),
        ('recovery_delay', 'float', 'The recovery delay (s) inserted between presynaptic pulses', {'index': True}),
        ('n_spikes_evoked', 'int', 'The number of presynaptic spikes evoked'),
    ],
    'test_pulse': [
        """A short, usually hyperpolarizing pulse used to test the resistance of pipette, cell access, or cell membrane.
        """,
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
        ('data', 'array', 'Numpy array of presynaptic recording sampled at '+_sample_rate_str),
        ('data_start_time', 'float', "Starting time of the data chunk, relative to the beginning of the recording"),
    ],
    'stim_spike': [
        "An action potential evoked by a stimulus pulse",
        ('pulse_id', 'stim_pulse.id', '', {'index': True}),
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
        ('data', 'array', 'numpy array of baseline data sampled at '+_sample_rate_str),
        ('mode', 'float', 'most common value in the baseline snippet'),
    ],
    'pulse_response': [
        "A chunk of postsynaptic recording taken during a presynaptic pulse stimulus",
        ('recording_id', 'recording.id', 'The full recording from which this pulse was extracted', {'index': True}),
        ('pulse_id', 'stim_pulse.id', 'The presynaptic pulse', {'index': True}),
        ('pair_id', 'pair.id', 'The pre-post cell pair involved in this pulse response', {'index': True}),
        ('start_time', 'float', 'Starting time of this chunk of the recording in seconds, relative to the beginning of the recording'),
        ('data', 'array', 'numpy array of response data sampled at '+_sample_rate_str),
        ('baseline_id', 'baseline.id'),
    ],
}




#----------- define ORM classes -------------

ORMBase = declarative_base()

class NDArray(TypeDecorator):
    """For marshalling arrays in/out of binary DB fields.
    """
    impl = LargeBinary
    
    def process_bind_param(self, value, dialect):
        buf = io.BytesIO()
        np.save(buf, value, allow_pickle=False)
        return buf.getvalue()
        
    def process_result_value(self, value, dialect):
        buf = io.BytesIO(value)
        return np.load(buf, allow_pickle=False)


class FloatType(TypeDecorator):
    """For marshalling float types (including numpy).
    """
    impl = Float
    
    def process_bind_param(self, value, dialect):
        if value is None:
            return None
        return float(value)
        
    #def process_result_value(self, value, dialect):
        #buf = io.BytesIO(value)
        #return np.load(buf, allow_pickle=False)


_coltypes = {
    'int': Integer,
    'float': FloatType,
    'bool': Boolean,
    'str': String,
    'date': Date,
    'datetime': DateTime,
    'array': NDArray,
    'object': JSONB,
}


def generate_mapping(table, schema):
    """Generate an ORM mapping class from an entry in table_schemas.
    """
    name = table.capitalize()
    table_args = {}
    if isinstance(schema[0], str):
        table_args['comment'] = schema[0]
        schema = schema[1:]
    
    props = {
        '__tablename__': table,
        '__table_args__': table_args,
        'id': Column(Integer, primary_key=True),
    }
    for column in schema:
        colname, coltype = column[:2]
        kwds = {} if len(column) < 4 else column[3]
        kwds['comment'] = None if len(column) < 3 else column[2]
            
        if coltype not in _coltypes:
            if not coltype.endswith('.id'):
                raise ValueError("Unrecognized column type %s" % coltype)
            props[colname] = Column(Integer, ForeignKey(coltype), **kwds)
        else:
            ctyp = _coltypes[coltype]
            props[colname] = Column(ctyp, **kwds)
    
    props['time_created'] = Column(DateTime, default=func.now())
    props['time_modified'] = Column(DateTime, onupdate=func.current_timestamp())
    props['meta'] = Column(JSONB)

    return type(name, (ORMBase,), props)


def _generate_mapping(table):
    return generate_mapping(table, table_schemas[table])


def create_all_mappings():
    global Slice, Experiment, Electrode, Cell, Pair, SyncRec, Recording, PatchClampRecording, MultiPatchProbe
    global TestPulse, StimPulse, StimSpike, PulseResponse, Baseline

    # Generate ORM mapping classes
    Slice = _generate_mapping('slice')
    Experiment = _generate_mapping('experiment')
    Electrode = _generate_mapping('electrode')
    Cell = _generate_mapping('cell')
    Pair = _generate_mapping('pair')
    SyncRec = _generate_mapping('sync_rec')
    Recording = _generate_mapping('recording')
    PatchClampRecording = _generate_mapping('patch_clamp_recording')
    MultiPatchProbe = _generate_mapping('multi_patch_probe')
    TestPulse = _generate_mapping('test_pulse')
    StimPulse = _generate_mapping('stim_pulse')
    StimSpike = _generate_mapping('stim_spike')
    PulseResponse = _generate_mapping('pulse_response')
    Baseline = _generate_mapping('baseline')

    # Set up relationships
    Slice.experiments = relationship("Experiment", order_by=Experiment.id, back_populates="slice")
    Experiment.slice = relationship("Slice", back_populates="experiments")

    Experiment.sync_recs = relationship(SyncRec, order_by=SyncRec.id, back_populates="experiment", cascade='delete', single_parent=True)
    SyncRec.experiment = relationship(Experiment, back_populates='sync_recs')

    Experiment.electrodes = relationship(Electrode, order_by=Electrode.id, back_populates="experiment", cascade="delete", single_parent=True)
    Electrode.experiment = relationship(Experiment, back_populates="electrodes")

    Electrode.cell = relationship(Cell, back_populates="electrode", cascade="delete", single_parent=True)
    Cell.electrode = relationship(Electrode, back_populates="cell", single_parent=True)

    Experiment.pairs = relationship(Pair, back_populates="experiment", cascade="delete", single_parent=True)
    Pair.experiment = relationship(Experiment, back_populates="pairs")

    Pair.pre_cell = relationship(Cell, foreign_keys=[Pair.pre_cell_id])
    #Cell.pre_pairs = relationship(Pair, back_populates="pre_cell", single_parent=True, foreign_keys=[Pair.pre_cell])

    Pair.post_cell = relationship(Cell, foreign_keys=[Pair.post_cell_id])
    #Cell.post_pairs = relationship(Pair, back_populates="post_cell", single_parent=True, foreign_keys=[Pair.post_cell])

    Electrode.recordings = relationship(Recording, back_populates="electrode", cascade="delete", single_parent=True)
    Recording.electrode = relationship(Electrode, back_populates="recordings")

    SyncRec.recordings = relationship(Recording, order_by=Recording.id, back_populates="sync_rec", cascade="delete", single_parent=True)
    Recording.sync_rec = relationship(SyncRec, back_populates="recordings")

    Recording.patch_clamp_recording = relationship(PatchClampRecording, back_populates="recording", cascade="delete", single_parent=True)
    PatchClampRecording.recording = relationship(Recording, back_populates="patch_clamp_recording")

    PatchClampRecording.multi_patch_probe = relationship(MultiPatchProbe, back_populates="patch_clamp_recording", cascade="delete", single_parent=True)
    MultiPatchProbe.patch_clamp_recording = relationship(PatchClampRecording, back_populates="multi_patch_probe")

    PatchClampRecording.nearest_test_pulse = relationship(TestPulse, cascade="delete", single_parent=True, foreign_keys=[PatchClampRecording.nearest_test_pulse_id])
    #TestPulse.patch_clamp_recording = relationship(PatchClampRecording)

    Recording.stim_pulses = relationship(StimPulse, back_populates="recording", cascade="delete", single_parent=True)
    StimPulse.recording = relationship(Recording, back_populates="stim_pulses")

    StimSpike.pulse = relationship(StimPulse, back_populates="spikes")
    StimPulse.spikes = relationship(StimSpike, back_populates="pulse", single_parent=True)

    Recording.baselines = relationship(Baseline, back_populates="recording", cascade="delete", single_parent=True)
    Baseline.recording = relationship(Recording, back_populates="baselines")

    PulseResponse.recording = relationship(Recording)
    PulseResponse.stim_pulse = relationship(StimPulse)
    Pair.pulse_responses = relationship(PulseResponse, back_populates='pair', single_parent=True)
    PulseResponse.pair = relationship(Pair, back_populates='pulse_responses')
    PulseResponse.baseline = relationship(Baseline)


#-------------- initial DB access ----------------

engine = create_engine(config.synphys_db_host + '/' + config.synphys_db)
# external users should create sessions from here.
Session = sessionmaker(bind=engine)

create_all_mappings()



def reset_db():
    """Drop the existing synphys database and initialize a new one.
    """
    pg_engine = create_engine(config.synphys_db_host + '/postgres')
    with pg_engine.begin() as conn:
        conn.connection.set_isolation_level(0)
        try:
            conn.execute('drop database %s' % config.synphys_db)
        except sqlalchemy.exc.ProgrammingError as err:
            if 'does not exist' not in err.message:
                raise

        conn.execute('create database %s' % config.synphys_db)

    # reconnect to DB
    global engine
    engine = create_engine(config.synphys_db_host + '/' + config.synphys_db)

    # Grant readonly permissions
    ro_user = config.synphys_db_readonly_user
    if ro_user is not None:
        with engine.begin() as conn:
            conn.execute('ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT ON TABLES TO %s;' % ro_user)
            # should only be needed if there are already tables present
            #conn.execute('GRANT SELECT ON ALL TABLES IN SCHEMA public TO %s;' % ro_user)
    
    # Create all tables
    global ORMBase
    ORMBase = declarative_base()
    create_all_mappings()
    ORMBase.metadata.create_all(engine)


def vacuum():
    with engine.begin() as conn:
        conn.connection.set_isolation_level(0)
        conn.execute('vacuum analyze')


def default_session(fn):
    def wrap_with_session(*args, **kwds):
        close = False
        if kwds.get('session', None) is None:
            kwds['session'] = Session()
            close = True
        try:
            ret = fn(*args, **kwds)
            return ret
        finally:
            if close:
                kwds['session'].close()
    return wrap_with_session    


@default_session
def slice_from_timestamp(ts, session=None):
    slices = session.query(Slice).filter(Slice.acq_timestamp==ts).all()
    if len(slices) == 0:
        raise KeyError("No slice found for timestamp %s" % ts)
    elif len(slices) > 1:
        raise KeyError("Multiple slices found for timestamp %s" % ts)
    
    return slices[0]


@default_session
def experiment_from_timestamp(ts, session=None):
    expts = session.query(Experiment).filter(Experiment.acq_timestamp==ts).all()
    if len(expts) == 0:
        raise KeyError("No experiment found for timestamp %s" % ts)
    elif len(expts) > 1:
        raise RuntimeError("Multiple experiments found for timestamp %s" % ts)
    
    return expts[0]



if __name__ == '__main__':
    # start a session
    session = Session()
    
    sl = Slice(lims_specimen_name="xxxxx", surface='medial')
    exp1 = Experiment(slice=sl, acsf='MP ACSF 1')
    exp2 = Experiment(slice=sl, acsf='MP ACSF 1')
    
    srec1 = SyncRec(experiment=exp1)
    srec2 = SyncRec(experiment=exp2)
    srec3 = SyncRec(experiment=exp2)
    
    rec1 = Recording(sync_rec=srec1)
    rec2 = Recording(sync_rec=srec2)
    rec3 = Recording(sync_rec=srec3)
    rec4 = Recording(sync_rec=srec3)
    
    pcrec1 = PatchClampRecording(recording=rec1)
    pcrec2 = PatchClampRecording(recording=rec2)
    pcrec3 = PatchClampRecording(recording=rec3)
    
    tp1 = TestPulse()
    tp2 = TestPulse()
    pcrec1.nearest_test_pulse = tp1
    pcrec2.nearest_test_pulse = tp2
    
    session.add_all([sl, exp1, exp2, srec1, srec2, srec3, rec1, rec2, rec3, rec4, pcrec1, pcrec2, pcrec3, tp1, tp2])
    session.commit()
