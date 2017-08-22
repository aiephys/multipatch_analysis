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
        ('surface', 'str', 'The surface of the slice facing up during the experiment (eg "medial")'),
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







def load_data(self, expt, pre=None, post=None):
    """Populate the database from raw data
    """
    expt_rows = self.tables['experiment']
    expt_index = self.tables['_expt_index']
    cell_rows = self.tables['cell']
    srec_rows = self.tables['sync_rec']
    rec_rows = self.tables['recording']
    tp_rows = self.tables['test_pulse']
    pulse_rows = self.tables['stim_pulse']
    spike_rows = self.tables['stim_spike']
    response_rows = self.tables['response']
    baseline_rows = self.tables['baseline']
    
    prof = Profiler(disabled=True, delayed=False)
    
    if expt.expt_id in expt_index:
        print("Cached: %s" % expt)
        raise NotImplementedError()
    
    prof.mark('start')
    expt_id = len(expt_rows)
    expt_rows.append({'id': expt_id, 'expt_key': expt.expt_id, 'internal_id': -1,
        'acsf_id': -1, 'temperature': np.nan, 'age': expt.age, 'genotype': None,
        'date': expt.date})
    expt_index[expt.expt_id] = expt_id
    
    cell_ids = {}
    for cell in expt.cells.values():
        cell_id = len(cell_rows)
        # mapping from experiment's internal ID for this cell to global cell ID 
        cell_ids[cell.cell_id] = cell_id
        cell_rows.append({'id': cell_id, 'expt_id': expt_id, 
            'device_key': cell.cell_id, 'cre_type': cell.cre_type,
            'pass_qc': cell.pass_qc, 'position': cell.position,
            'depth': cell.depth})
    prof.mark('cells')


    expt_data = expt.data
    for srec in expt_data.contents:
        srec_id = len(srec_rows)
        srec_rows.append({'id': srec_id, 'expt_id': expt_id,
            'sync_rec_key': srec.key})
        rec_key_id_map = {}
        pulse_key_n_id_map = {}
        for rec in srec.recordings:
            rec_id = len(rec_rows)
            rec_key_id_map[rec.device_id] = rec_id
            tp_id = len(tp_rows)
            cell_id = cell_ids[rec.device_id + 1]
            psa = PulseStimAnalyzer.get(rec)
            ind_freq, recovery_delay = psa.stim_params()
            
            rec_rows.append({'id': rec_id, 'sync_rec_id': srec_id, 'cell_id': cell_id,
                'device_key': rec.device_id, 'stim_name': rec.meta['stim_name'],
                'clamp_mode': rec.clamp_mode, 'test_pulse_id': tp_id,
                'sample_rate': rec['primary'].sample_rate,
                'induction_freq': ind_freq, 'recovery_delay': recovery_delay})
            
            pulses = psa.pulses()
            if pulses[0][2] < 0:
                # test pulse
                tp = pulses[0]
            else:
                tp = (None, None)
            tp_rows.append({'id': tp_id, 'recording_id': rec_id,
                'pulse_start': tp[0], 'pulse_stop': tp[1],
                })
            
            if pre is None or rec.device_id == pre:
                pulse_n_id_map = {}
                for i,pulse in enumerate(pulses):
                    pulse_id = len(pulse_rows)
                    pulse_n_id_map[i] = pulse_id
                    pulse_key_n_id_map[(rec.device_id, i)] = pulse_id
                    pulse_rows.append({'id': pulse_id, 'recording_id': rec_id,
                        'pulse_number': i, 'onset_index': pulse[0],
                        'length': pulse[1]-pulse[0], 'amplitude': pulse[2],
                        'n_spikes': 0})
            
                spikes = psa.evoked_spikes()
                for sp in spikes:
                    sp_id = len(spike_rows)
                    pulse_id = pulse_n_id_map[sp['pulse_n']]
                    srow = {'id': sp_id, 'recording_id': rec_id, 'pulse_id': pulse_id}
                    pulse_rows[pulse_id]['n_spikes'] += 1
                    if sp['spike'] is not None:
                        srow.update(sp['spike'])
                    spike_rows.append(srow)
            
        mpa = MultiPatchSyncRecAnalyzer(srec)
        for pre_dev in srec.devices:
            if pre is not None and pre_dev != pre:
                continue
            
            for post_dev in srec.devices:
                if post is not None and post_dev != post:
                    continue
                
                responses = mpa.get_spike_responses(srec[pre_dev], srec[post_dev], align_to='pulse', require_spike=False)
                for resp in responses:
                    resp_id = len(response_rows)
                    bl_id = len(baseline_rows)
                    baseline_rows.append({'id': bl_id,
                        'recording_id': rec_key_id_map[post_dev],
                        'start_index': resp['baseline_start'],
                        'stop_index': resp['baseline_stop'],
                        'data': resp['baseline'].downsample(f=50000).data,
                        'value': float_mode(resp['baseline'].data),
                    })
                    response_rows.append({'id': resp_id, 
                        'recording_id': rec_key_id_map[post_dev],
                        'pulse_id': pulse_key_n_id_map[(pre_dev, resp['pulse_n'])],
                        'start_index': resp['rec_start'], 'stop_index': resp['rec_stop'],
                        'baseline_id': bl_id,
                        'data': resp['response'].downsample(f=50000).data,
                    })
    

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
