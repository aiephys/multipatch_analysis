"""
Accumulate all experiment data into a set of linked tables.
"""
import os, pickle, resource
import numpy as np
from pandas import DataFrame

from neuroanalysis.baseline import float_mode
from connection_detection import PulseStimAnalyzer, MultiPatchSyncRecAnalyzer


class ExperimentDatabase(object):
    
    table_schemas = {
        'experiment': [
            ('expt_key', 'object'),
            ('internal_id', 'int'),
            ('acsf_id', 'int'),
            ('temperature', 'float'),
            ('age', 'int'),
            ('genotype', 'string'),
            ('date', 'object')
        ],
        'cell': [
            ('expt_id', 'int'),
            ('cre_type', 'string'),
            ('device_key', 'object'),
            ('patch_start', 'float'),
            ('patch_stop', 'float'),
            ('initial_seal_resistance', 'float'),
            ('initial_pipette_resistance', 'float'),
            ('final_pipette_resistance', 'float'),
            ('has_biocytin', 'bool'),
            ('has_dye_fill', 'bool'),
            ('pass_qc', 'bool'),
            ('pass_spike_qc', 'bool'),
        ],
        'sync_rec': [
            ('expt_id', 'int'),
            ('sync_rec_key', 'object'),
            ('time_post_patch', 'float'),
        ],
        'recording': [
            ('sync_rec_id', 'int'),
            ('device_key', 'object'),
            ('stim_name', 'object'),
            ('induction_freq', 'float'),
            ('recovery_delay', 'float'),
            ('clamp_mode', 'string'),
            ('test_pulse_id', 'int'),
            ('sample_rate', 'float'),
        ],
        'test_pulse': [
            ('recording_id', 'int'),
            ('pulse_start', 'int'),
            ('pulse_stop', 'int'),
            ('baseline_current', 'float'),
            ('baseline_voltage', 'float'),
            ('input_resistance', 'float'),
            ('access_resistance', 'float'),
        ],
        'stim_pulse': [
            ('recording_id', 'int'),
            ('pulse_number', 'int'),
            ('onset_time', 'float'),
            ('onset_index', 'int'),
            ('amplitude', 'float'),
            ('length', 'int'),
            ('n_spikes', 'int'),                           # number of spikes evoked
        ],
        'stim_spike': [                                  # One evoked action potential
            ('recording_id', 'int'),
            ('pulse_id', 'int'),
            ('peak_index', 'int'),
            ('peak_diff', 'float'),
            ('rise_index', 'int'),
            ('max_dvdt', 'float'),
        ],
        'response': [                                    # One evoked synaptic response
            ('recording_id', 'int'),
            ('pulse_id', 'int'),
            ('start_index', 'int'),
            ('stop_index', 'int'),
            ('data', 'object'),
            ('baseline_id', 'int'),
        ],
        'baseline': [                                    # One snippet of baseline data
            ('recording_id', 'int'),
            ('start_index', 'int'),                        # start/stop indices of baseline snippet
            ('stop_index', 'int'),                         #   relative to full recording
            ('data', 'object'),                             # array containing baseline snippet
            ('value', 'float'),                             # median or mode baseline value
        ],
    }
        

    def __init__(self, cache=None):
        self.tables = {name:[] for name in self.table_schemas}
        self.tables['_expt_index'] = {}
        self.tables['_version'] = 1
        
        self.cache = cache
        if cache is not None and os.path.isfile(cache):
            self.load_cache()
            
    def load_cache(self):
        self.tables.update(pickle.load(open(self.cache, 'rb')))
        
    def store_cache(self):
        pickle.dump(self.tables, open(self.cache, 'wb'))

    def load_data(self, expts):
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
        
        for expt in expts:
            if expt.expt_id in expt_index:
                print("Cached: %s" % expt)
                continue
            
            print("Load: %s" % expt)
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
            expt_data = expt.data
            for srec in expt_data.contents:
                srec_id = len(srec_rows)
                srec_rows.append([{'id': srec_id, 'expt_id': expt_id,
                    'sync_rec_key': srec.key}])
                rec_key_id_map = {}
                pulse_key_n_id_map = {}
                for rec in srec.recordings:
                    rec_id = len(rec_rows)
                    rec_key_id_map[rec.device_id] = rec_id
                    tp_id = len(tp_rows)
                    rec_rows.append({'id': rec_id, 'sync_rec_id': srec_id,
                        'device_key': rec.device_id, 'stim_name': rec.meta['stim_name'],
                        'clamp_mode': rec.clamp_mode, 'test_pulse_id': tp_id,
                        'sample_rate': rec['primary'].sample_rate})
                    
                    psa = PulseStimAnalyzer.get(rec)
                    pulses = psa.pulses()
                    if pulses[0][2] < 0:
                        # test pulse
                        tp = pulses[0]
                    else:
                        tp = (None, None)
                    tp_rows.append({'id': tp_id, 'recording_id': rec_id,
                        'pulse_start': tp[0], 'pulse_stop': tp[1],
                        })
                        
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
                    for post_dev in srec.devices:
                        responses = mpa.get_spike_responses(srec[pre_dev], srec[post_dev], align_to='pulse', require_spike=False)
                        for resp in responses:
                            resp_id = len(response_rows)
                            bl_id = len(baseline_rows)
                            baseline_rows.append({'id': bl_id,
                                'recording_id': rec_key_id_map[post_dev],
                                'start_index': resp['baseline_start'],
                                'stop_index': resp['baseline_stop'],
                                'data': resp['baseline'].data,
                                'value': float_mode(resp['baseline'].data),
                            })
                            response_rows.append({'id': resp_id, 
                                'recording_id': rec_key_id_map[post_dev],
                                'pulse_id': pulse_key_n_id_map[(pre_dev, resp['pulse_n'])],
                                'start_index': resp['rec_start'], 'stop_index': resp['rec_stop'],
                                'baseline_id': bl_id, 'data': resp['response'].data,
                            })
                        
            expt.close_data()
            
            self.store_cache()
            
            mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            print("memory usage: %d" % mem)
            if mem > 16000000:
                print("   it is a good day to die.")
                sys.exit(-1)

        
    def make_table(self, name, data):
        schema = self.table_schemas[name]
        schema.insert(0, ('id', 'uint'))
        a = np.empty(0, schema)
        table = DataFrame(a)
        table = table.append(data)
        return table
            
    def get_table(self, name):
        if name not in self.tables:
            a = np.empty(0, self.table_schemas[name])
            self.tables[name] = DataFrame(a)
        return self.tables[name]
        

if __name__ == '__main__':
    import pyqtgraph as pg
    pg.dbg()
    
    from experiment_list import ExperimentList
    all_expts = ExperimentList(cache='expts_cache.pkl')
    
    
    db_cache_file = 'database.pkl'
    db = ExperimentDatabase(cache=db_cache_file)
    db.load_data(all_expts)
    db.store_cache()


