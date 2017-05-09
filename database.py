"""
Accumulate all experiment data into a set of linked tables.


Response Table : one row per stim/response

- Response info:
    link to data
    response snippet (to allow mass averaging)
    psp fit params
- Baseline info:
    link to data
    snippet
- Test pulse info:
    link to data
    holding current
    baseline potential
    input resistance
    access resistance
- Stim info
    holding
    clamp mode
    ind freq
    rec delay
- Expt info
    expt_id
    sweep_id
    cell_ids
    cre types
    connected
    time post-patching
    
"""
import os, pickle
import numpy as np
from pandas import DataFrame


class ExperimentDatabase(object):
    
    table_schemas = {
        'experiment': [
            ('expt_key', 'object'),
            ('internal_id', 'int'),
            ('acsf_id', 'int'),
            ('temperature', 'float'),
            ('age', 'int'),
            ('genotype', 'string'),
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
            ('peak_latency', 'float'),
            ('peak_amplitude', 'float'),
            ('onset_latency', 'float'),
            ('onset_amplitude', 'float'),
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
        self.tables = {}
        self.cache = cache
        if cache is not None and os.path.isfile(cache):
            self.load_cache()
            
    def load_cache(self):
        self.data = pickle.load(open(self.cache, 'rb'))
        
    def store_cache(self):
        pickle.dump(self.tables, open(self.cache, 'wb'))

    def load_data(self, expts):
        """Populate the database from raw data
        """
        expt_recs = []
        for expt in expts:
            expt_recs.append({'expt_key': expt.expt_id, 'internal_id': -1,
                              'acsf_id': -1, 'temperature': np.nan, 'age': expt.age,
                              'genotype': None})
            
        self.tables['experiment'] = self.make_table('experiment', expt_recs)
        
    def make_table(self, name, data):
        a = np.empty(0, self.table_schemas[name])
        table = DataFrame(a)
        return table.append(data)
            
    def get_table(self, name):
        if name not in self.tables:
            a = np.empty(0, self.table_schemas[name])
            self.tables[name] = DataFrame(a)
        return self.tables[name]
        

if __name__ == '__main__':
    from experiment_list import ExperimentList
    all_expts = ExperimentList(cache='expts_cache.pkl')
    
    
    db_cache_file = 'database.pkl'
    db = ExperimentDatabase(cache=db_cache_file)
    db.load_data(all_expts)
    db.store_cache()


