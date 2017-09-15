"""
Accumulate all experiment data into a set of linked tables.
"""
import io
import numpy as np

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Boolean, Float, Date, DateTime, LargeBinary, ForeignKey
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship
from sqlalchemy.types import TypeDecorator
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.expression import func


from neuroanalysis.baseline import float_mode
from connection_detection import PulseStimAnalyzer, MultiPatchSyncRecAnalyzer
from config import synphys_db


table_schemas = {
    'slice': [
        "All brain slices on which an experiment was attempted.",
        ('acq_timestamp', 'datetime', 'Creation timestamp for slice data acquisition folder.'),
        ('species', 'str', 'Human | mouse (from LIMS)'),
        ('age', 'int', 'Specimen age (in days) at time of dissection (from LIMS)'),
        ('genotype', 'str', 'Specimen donor genotype (from LIMS)'),
        ('orientation', 'str', 'Orientation of the slice plane (eg "sagittal"; from LIMS specimen name)'),
        ('surface', 'str', 'The surface of the slice exposed during the experiment (eg "left"; from LIMS specimen name)'),
        ('hemisphere', 'str', 'The brain hemisphere from which the slice originated. (from LIMS specimen name)'),
        ('quality', 'int', 'Experimenter subjective slice quality assessment (0-5)'),
        ('slice_time', 'datetime', 'Time when this specimen was sliced'),
        ('slice_conditions', 'object', 'JSON containing solutions, perfusion, incubation time, etc.'),
        ('lims_specimen_name', 'str', 'Name of LIMS "slice" specimen'),
        ('original_path', 'str', 'Original path of the slice folder on the acquisition rig'),
        ('submission_data', 'object'),          # structure generated for original submission
    ],
    'experiment': [
        "A group of cells patched simultaneously in the same slice.",
        ('original_path', 'str', 'Describes original location of raw data'),
        ('acq_timestamp', 'datetime', 'Creation timestamp for site data acquisition folder.'),
        ('slice_id', 'slice.id'),
        ('target_region', 'str', 'The intended brain region for this experiment'),
        ('internal', 'str', 'The name of the internal solution used in this experiment. '
                            'The solution should be described in the pycsf database.'),
        ('acsf', 'str', 'The name of the ACSF solution used in this experiment. '
                        'The solution should be described in the pycsf database.'),
        ('target_temperature', 'float'),
        ('date', 'datetime'),
        ('lims_specimen_id', 'int', 'ID of LIMS "CellCluster" specimen.'),
        ('submission_data', 'object', 'structure generated for original submission.'),
        ('lims_trigger_id', 'int', 'ID used to query status of LIMS upload.'),
        ('connectivity_analysis_complete', 'bool'),
        ('kinetics_analysis_complete', 'bool'),
    ],
    'electrode': [
        "Each electrode records a patch attempt, whether or not it resulted in a "
        "successful cell recording.",
        ('expt_id', 'experiment.id'),
        ('patch_status', 'str', 'no seal, low seal, GOhm seal, tech fail, ...'),
        ('device_key', 'int'),
        ('initial_resistance', 'float'),
        ('initial_current', 'float'),
        ('pipette_offset', 'float'),
        ('final_resistance', 'float'),
        ('final_current', 'float'),
        ('notes', 'str'),
    ],
    'cell': [
        ('electrode_id', 'electrode.id'),
        ('cre_type', 'str'),
        ('patch_start', 'float'),
        ('patch_stop', 'float'),
        ('seal_resistance', 'float'),
        ('has_biocytin', 'bool'),
        ('has_dye_fill', 'bool'),
        ('pass_qc', 'bool'),
        ('pass_spike_qc', 'bool'),
        ('depth', 'float'),
        ('position', 'object'),
    ],
    
    'pair': [
        "All possible putative synaptic connections",
        ('pre_cell', 'cell.id'),
        ('post_cell', 'cell.id'),
        ('synapse', 'bool', 'Whether the experimenter thinks there is a synapse'),
        ('electrical', 'bool', 'whether the experimenter thinks there is a gap junction'),
    ],
    'sync_rec': [
        ('expt_id', 'experiment.id'),
        ('sync_rec_key', 'object'),
        ('meta', 'object'),
    ],
    'recording': [
        ('sync_rec_id', 'sync_rec.id', 'References the synchronous recording to which this recording belongs.'),
        ('device_key', 'object', 'Identifies the device that generated this recording (this is usually the MIES AD channel)'),
        ('start_time', 'datetime', 'The clock time at the start of this recording'),
    ],
    'patch_clamp_recording': [
        "Extra data for recordings made with a patch clamp amplifier",
        ('recording_id', 'recording.id'),
        ('electrode_id', 'electrode.id', 'References the patch electrode that was used during this recording'),
        ('clamp_mode', 'str', 'The mode used by the patch clamp amplifier: "ic" or "vc"'),
        ('patch_mode', 'str', "The state of the membrane patch. E.g. 'whole cell', 'cell attached', 'loose seal', 'bath', 'inside out', 'outside out'"),
        ('stim_name', 'object', "The name of the stimulus protocol"),
        ('baseline_potential', 'float'),
        ('baseline_current', 'float'),
        ('baseline_rms_noise', 'float'),
        #('nearest_test_pulse_id', 'test_pulse.id', "Reference to the test pulse for this recording"),
        #('lowpass_test_pulse_id', 'test_pulse.id', "Reference to the low-passed test pulse values"),
    ],
    'test_pulse': [
        ('patch_clamp_recording_id', 'patch_clamp_recording.id'),
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
        ('recording_id', 'recording.id'),
        ('pulse_number', 'int'),
        ('onset_time', 'float'),
        ('onset_index', 'int'),
        ('next_pulse_index', 'int'),      # index of the next pulse on any channel in the sync rec
        ('amplitude', 'float'),
        ('length', 'int'),
        ('n_spikes', 'int'),                           # number of spikes evoked
    ],
    'stim_spike': [
        "An evoked action potential",
        ('recording_id', 'recording.id'),
        ('pulse_id', 'stim_pulse.id'),
        ('peak_index', 'int'),
        ('peak_diff', 'float'),
        ('peak_val', 'float'),
        ('rise_index', 'int'),
        ('max_dvdt', 'float'),
    ],
    'baseline': [
        "A snippet of baseline data, matched to a postsynaptic recording",
        ('recording_id', 'recording.id', 'The recording from which this baseline snippet was extracted.'),
        ('start_index', 'int', 'start index of this snippet, relative to the beginning of the recording'),
        ('stop_index', 'int', 'stop index of this snippet, relative to the beginning of the recording'),
        ('data', 'array', 'array containing the baseline snippet'),
        ('mode', 'float', 'most common value in the baseline snippet'),
    ],
    'pulse_response': [
        "A postsynaptic recording taken during a presynaptic stimulus",
        ('recording_id', 'recording.id'),
        ('pulse_id', 'stim_pulse.id'),
        ('pair_id', 'pair.id'),
        ('start_index', 'int'),
        ('stop_index', 'int'),
        ('data', 'array'),
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


def _generate_mapping(table):
    """Generate an ORM mapping class from an entry in table_schemas.
    """
    name = table.capitalize()
    schema = table_schemas[table]
    table_args = {}
    if isinstance(schema[0], str):
        table_args['comment'] = schema[0]
        schema = schema[1:]
    
    props = {
        '__tablename__': table,
        '__table_args__': table_args,
        'id': Column(Integer, primary_key=True),
        'time_created': Column(DateTime, default=func.now()),
        'time_modified': Column(DateTime, onupdate=func.current_timestamp()),
    }
    for column in schema:
        colname, coltype = column[:2]
        if len(column) > 2:
            comment = column[2]
        else:
            comment = None
            
        if coltype not in _coltypes:
            if not coltype.endswith('.id'):
                raise ValueError("Unrecognized column type %s" % coltype)
            props[colname] = Column(Integer, ForeignKey(coltype), comment=comment)
        else:
            ctyp = _coltypes[coltype]
            props[colname] = Column(ctyp, comment=comment)
    return type(name, (ORMBase,), props)


# Generate ORM mapping classes
Slice = _generate_mapping('slice')
Experiment = _generate_mapping('experiment')
Electrode = _generate_mapping('electrode')
Cell = _generate_mapping('cell')
Pair = _generate_mapping('pair')
SyncRec = _generate_mapping('sync_rec')
Recording = _generate_mapping('recording')
PatchClampRecording = _generate_mapping('patch_clamp_recording')
TestPulse = _generate_mapping('test_pulse')
StimPulse = _generate_mapping('stim_pulse')
StimSpike = _generate_mapping('stim_spike')
PulseResponse = _generate_mapping('pulse_response')
Baseline = _generate_mapping('baseline')

# Set up relationships
Experiment.slice = relationship("Slice", back_populates="experiments")
Slice.experiments = relationship("Experiment", order_by=Experiment.id, back_populates="slice")

SyncRec.experiment = relationship(Experiment)
#Experiment.sync_recs = relationship("SyncRec", order_by=SyncRec.id, back_populates="experiment")

Recording.sync_rec = relationship(SyncRec)
PatchClampRecording.recording = relationship(Recording)
TestPulse.patch_clamp_recording = relationship(PatchClampRecording)

StimPulse.recording = relationship(Recording)

StimSpike.recording = relationship(Recording)
StimSpike.pulse = relationship(StimPulse)

Baseline.recording = relationship(Recording)

PulseResponse.recording = relationship(Recording)
PulseResponse.stim_pulse = relationship(StimPulse)
PulseResponse.baseline = relationship(Baseline)


#-------------- initial DB access ----------------

# connect to DB
engine = create_engine(synphys_db)


# recreate all tables in DB
# (just for initial development)
import sys
if '--reset-db' in sys.argv:
    ORMBase.metadata.drop_all(engine)
    ORMBase.metadata.create_all(engine)


# external users should create sessions from here.
Session = sessionmaker(bind=engine)







    

def slice_from_timestamp(ts):
    session = Session()
    slices = session.query(Slice).filter(Slice.acq_timestamp==ts).all()
    if len(slices) == 0:
        raise KeyError("No slice found for timestamp %s" % ts)
    elif len(slices) > 1:
        raise KeyError("Multiple slices found for timestamp %s" % ts)
    
    return slices[0]



if __name__ == '__main__':
    # start a session
    session = Session()
    
    sl = Slice(lims_specimen_name="xxxxx", surface='medial')
    exp = Experiment(slice=sl, acsf='MP ACSF 1')
    exp2 = Experiment(slice=sl, acsf='MP ACSF 1')
    
    session.add(sl)
    session.commit()
