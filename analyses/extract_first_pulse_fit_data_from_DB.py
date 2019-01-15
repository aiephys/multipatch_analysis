# Extract relevant data from the database for analysis
# NOTE: that you need to make sure the import on line 10,
# the query on line 35, the image path on line 68, and
# the out put file path on line 72 are appropriately 
# matched

from multipatch_analysis.database import database as db
import pandas as pd
import find_image_file
from multipatch_analysis.first_pulse_fits_average import AFPFForceSignJitterLatency   

# Initialize output dictionary
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

# Do the query. Note this can take several minutes
session=db.Session()
data=session.query(AFPFForceSignJitterLatency, db.Pair).join(db.Pair).all() #this need to correspond to import

# extract the relevant data from the query
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

# put data in dataframe
df=pd.DataFrame(data_dict)
df['uid']=df['uid'].astype(str)

# add the image path
df['image_path']=df.apply(lambda row: 
                          find_image_file.find_image_file('/home/corinnet/workspace/DBfit_pics/dynamic_weight_jitter_latency2018-11-12', #make sure this is correct
                          row.uid, str(row.pre_cell_id), str(row.post_cell_id)), axis=1)

# save data to csv                          
#df.to_csv('dynamic_weight_jitter_latency2018-11-12.csv')  #comment out after using so dont overwrite something