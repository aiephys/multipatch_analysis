# this code will enable querying from the *avg_first_pulse_fit_dynamic_w_latency_jitter2* 
# table in synphys9.  It will query the *avg_first_pulse_fit_dynamic_w_latency_jitter2* 
# table along with the *pair* table and save desired information into a csv file.

from multipatch_analysis.database import database as db
from multipatch_analysis.database import database as db
#import multipatch_analysis.connection_strength as cs 
from multipatch_analysis.database.database import TableGroup
import pandas as pd
#----------------------------------------------------------------
#------------Schema for accessing table--------------------------
#----------------------------------------------------------------

class FirstPulseFitTableGroup(TableGroup):
    """Fits first pulse for each individual sweeps.
    """
    schemas = {
        'avg_first_pulse_fit_dynamic_w_latency_jitter2': [
             """Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed qc in current clamp. During the fit the latency is forced
            to be the value +/-.5 ms found via fitting all of the pulses (available in the 
            connection_strength.ic_fit_xoffset). The heavily weighted section meant to 
            place more importance of the wave form during the rise time is shifted to 
            begin at the latency. Created via first_pulse_fits_average.py. 
            All units in SI.""",
            ('pair_id', 'pair.id', '', {'index': True}),
            ('uid', 'float','timestamp attached to the experiment for ease of viewing'),
            ('pre_cell_id', 'int', 'the integer id of the pre synaptic cell (from cell table)'),
            ('post_cell_id', 'int', 'the integer id of the post synaptic cell (from cell table)'),
            ('amp', 'float', 'amplitude '),
            ('latency', 'float', 'time elapsed since the time of presynaptic spike (max dv/dt)'),
            ('rise_time', 'float', 'rise time of psp', ),
            ('decay_tau', 'float', 'decay of psp'),
            ('avg_psp', 'array', 'array of the best fit voltage waveform starting 10 ms before pre-synaptic spike'),
            ('dt', 'float', 'time step of *avg_psp* array'),
            ('n_sweeps', 'int', 'number of sweeps used in the fit'),
            ('pulse_ids', 'object', 'data base pulse ids included in fit'),
            ('distance', 'float', 'distance between pairs'),
            ('NRMSE', 'float', 'error of fit'),
            ('synapse_sign', 'str', '"ex" or "in" (also in connection strength table but here for convenience)'),
            ('measured_baseline', 'float', 'average voltage measured between 10 and 1 ms before a spike'),
            ('measured_amp', 'float', 'amplitude within a window of 0.5 ms after spike initiation (max dv/dt) until end of array specified in the pulse_response table'),
            ('connected', 'bool', 'specifies whether human thinks their is a connection'),
            ],
        }
    
    def create_mappings(self):
        TableGroup.create_mappings(self)
        AFPFForceSignJitterLatency = self['avg_first_pulse_fit_dynamic_w_latency_jitter2']
        db.Pair.avg_first_pulse_fit_dynamic_w_latency_jitter2 = db.relationship(AFPFForceSignJitterLatency, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        AFPFForceSignJitterLatency.pair = db.relationship(db.Pair, back_populates="avg_first_pulse_fit_dynamic_w_latency_jitter2", single_parent=True)

first_pulse_fit_tables = FirstPulseFitTableGroup()

def init_tables():
    global AFPFForceSignJitterLatency
    AFPFForceSignJitterLatency = first_pulse_fit_tables['avg_first_pulse_fit_dynamic_w_latency_jitter2']

init_tables()

#----------------------------------------------------------------
#------------query and put in a dataframe--------------------------
#----------------------------------------------------------------
'''Note that this query can take several minutes'''

# initialize dictionary
data_dict= {'uid':[],
            'pre_cell_id':[],
            'post_cell_id':[],
            'pre_cre':[],
            'post_cre':[],
            'amp':[], 
            'NRMSE':[], 
            'decay_tau':[], 
            'latency':[], 
            'rise_time':[], 
            'syn_excitation':[],
            'distance':[],
            'boolean_connection':[],
            'acsf':[],
            'measured_amp':[],
            'measured_baseline':[],
            'n_sweeps':[],
            'pre_layer':[],
            'post_layer':[]}

# do query
session=db.Session()
data=session.query(AFPFForceSignJitterLatency, db.Pair).join(db.Pair).all() #this need to correspond to import

# extract the relevant data from the query and place in dictionary
for afpf, pair in data:
    #stuff from average_first_pulse_fit table
    data_dict['uid'].append(afpf.__dict__['uid'])
    data_dict['pre_cell_id'].append(afpf.__dict__['pre_cell_id'])
    data_dict['post_cell_id'].append(afpf.__dict__['post_cell_id'])
    data_dict['amp'].append(afpf.__dict__['amp'])
    data_dict['NRMSE'].append(afpf.__dict__['NRMSE'])
    data_dict['decay_tau'].append(afpf.__dict__['decay_tau'])
    data_dict['rise_time'].append(afpf.__dict__['rise_time'])
    data_dict['latency'].append(afpf.__dict__['latency'])
    data_dict['boolean_connection'].append(afpf.__dict__['connected'])
    data_dict['syn_excitation'].append(afpf.__dict__['synapse_sign'])
    data_dict['measured_amp'].append(afpf.__dict__['measured_amp'])
    data_dict['measured_baseline'].append(afpf.__dict__['measured_baseline'])
    data_dict['n_sweeps'].append(afpf.__dict__['n_sweeps'])
    #stuff from pair table
    data_dict['pre_cre'].append(pair.pre_cell.cre_type)
    data_dict['post_cre'].append(pair.post_cell.cre_type)
    data_dict['distance'].append(pair.distance)
    data_dict['acsf'].append(pair.experiment.acsf)
    data_dict['pre_layer'].append(pair.pre_cell.target_layer)
    data_dict['post_layer'].append(pair.pre_cell.target_layer)
session.close()

#put dictionary into a dataframe
df=pd.DataFrame(data_dict)
df['uid']=df['uid'].astype(str)

#save data frame
filename='my_sweet_file.csv'
df.to_csv(filename)  #comment out after using so dont overwrite something