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
    'slice': [                            # most of this should be pulled from external sources
        ('acq_timestamp', 'datetime', 'Creation timestamp for slice data acquisition folder.'),
        ('species', 'str', 'Human | mouse'),
        ('age', 'int', 'Specimen age (in days) at time of dissection'),
        ('genotype', 'str', 'Specimen donor genotype (from LIMS)'),
        ('orientation', 'str', 'Orientation of the slice plane (eg "sagittal")'),
        ('surface', 'str', 'The surface of the slice exposed during the experiment (eg "left")'),
        ('hemisphere', 'str', 'The brain hemisphere from which the slice originated.'),
        ('quality', 'int', 'Experimenter subjective slice quality assessment (0-5)'),
        ('slice_time', 'datetime', 'Time when this specimen was sliced'),
        ('slice_conditions', 'object', 'JSON containing solutions, perfusion, incubation time, etc.'),
        ('lims_specimen_name', 'str', 'Name of LIMS "slice" specimen'),
        ('original_path', 'str', 'Original path of the slice folder on the acquisition rig'),
        ('submission_data', 'object'),          # structure generated for original submission
    ],
    'experiment': [
        "A group of cells patched simultaneously in the same slice.",
        ('original_path', 'object', 'Describes original location of raw data'),
        ('slice_id', 'slice.id'),
        ('region', 'str', 'The intended brain region for this experiment'),
        ('internal', 'str', 'The name of the internal solution used in this experiment. '
                            'The solution should be described in the pycsf database.'),
        ('acsf', 'str', 'The name of the ACSF solution used in this experiment. '
                        'The solution should be described in the pycsf database.'),
        ('temperature', 'float'),
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
        ('patch_status', 'str'),
        ('device_key', 'int'),
    ],
    'cell': [
        ('electrode_id', 'electrode.id'),
        ('cre_type', 'str'),
        ('patch_start', 'float'),
        ('patch_stop', 'float'),
        ('initial_seal_resistance', 'float'),
        ('initial_pipette_resistance', 'float'),
        ('final_pipette_resistance', 'float'),
        ('has_biocytin', 'bool'),
        ('has_dye_fill', 'bool'),
        ('pass_qc', 'bool'),
        ('pass_spike_qc', 'bool'),
        ('depth', 'float'),
        ('position', 'object'),
    ],
    
    'pair': [     # table of all POSSIBLE connections
        ('pre_cell', 'cell.id'),
        ('post_cell', 'cell.id'),
        ('synapse', 'bool'),       # Whether the experimenter thinks there is a synapse
        ('electrical', 'bool'),    # whether the experimenter thinks there is a gap junction
    ],
        # NOTE: add individual per-pair analyses to new tables.
    
    'sync_rec': [
        ('expt_id', 'experiment.id'),
        ('sync_rec_key', 'object'),
        ('time_post_patch', 'float'),
    ],
    'recording': [
        ('sync_rec_id', 'sync_rec.id'),
        ('device_key', 'object'),
        ('electrode_id', 'electrode.id'),
        ('clamp_mode', 'object'),
        ('stimulus', 'object'),   # contains name, induction freq, recovery delay
        ('test_pulse', 'object'),  # contains pulse_start, pulse_stop, baseline_current,
                                   # baseline_voltage, input_resistance, access_resistance
        ('sample_rate', 'float'),
    ],
    'stim_pulse': [
        ('recording_id', 'recording.id'),
        ('pulse_number', 'int'),
        ('onset_time', 'float'),
        ('onset_index', 'int'),
        ('next_pulse_index', 'int'),      # index of the next pulse on any channel in the sync rec
        ('amplitude', 'float'),
        ('length', 'int'),
        ('n_spikes', 'int'),                           # number of spikes evoked
    ],
    'stim_spike': [                                  # One evoked action potential
        ('recording_id', 'recording.id'),
        ('pulse_id', 'stim_pulse.id'),
        ('peak_index', 'int'),
        ('peak_diff', 'float'),
        ('peak_val', 'float'),
        ('rise_index', 'int'),
        ('max_dvdt', 'float'),
    ],
    'pulse_response': [                                    # One evoked synaptic response
        ('recording_id', 'recording.id'),
        ('pulse_id', 'stim_pulse.id'),
        ('pair_id', 'pair.id'),
        ('start_index', 'int'),
        ('stop_index', 'int'),
        ('data', 'object'),
        ('baseline_id', 'baseline.id'),
    ],
    'baseline': [                                    # One snippet of baseline data
        ('recording_id', 'recording.id'),
        ('start_index', 'int'),                        # start/stop indices of baseline snippet
        ('stop_index', 'int'),                         #   relative to full recording
        ('data', 'bytes'),                             # array containing baseline snippet
        ('value', 'float'),                             # median or mode baseline value
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


_coltypes = {
    'int': Integer,
    'float': Float,
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

# Set up relationships
Experiment.slice = relationship("Slice", back_populates="experiments")
Slice.experiments = relationship("Experiment", order_by=Experiment.id, back_populates="slice")



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
    return slices



if __name__ == '__main__':
    # start a session
    session = Session()
    
    sl = Slice(lims_specimen_name="xxxxx", surface='medial')
    exp = Experiment(slice=sl, acsf='MP ACSF 1')
    exp2 = Experiment(slice=sl, acsf='MP ACSF 1')
    
    session.add(sl)
    session.commit()
