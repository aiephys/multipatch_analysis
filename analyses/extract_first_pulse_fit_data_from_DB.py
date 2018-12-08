from multipatch_analysis.database import database as db
#from multipatch_analysis.database.database import TableGroup
#import multipatch_analysis.connection_strength as cs 
import pandas as pd
import find_image_file

#THIS DOES QUERY, EXTRACTS DATA, AND SAVES IN A CSV.
# MAKE SURE YOU HAVE THE FORCE VERSUS ANY CORRECT THOUGHOUT THIS SECTION:
# search for any or force and make sure all are converted
'''Note that this query can take several minutes'''
from multipatch_analysis.first_pulse_fits_average import VClamp   
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

#Note that it is annoying that 
session=db.Session()
data=session.query(VClamp, db.Pair).join(db.Pair).all() #this need to correspond to import
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
df=pd.DataFrame(data_dict)
df['uid']=df['uid'].astype(str)
df['image_path']=df.apply(lambda row: 
                          find_image_file.find_image_file('/home/corinnet/workspace/DBfit_pics/vclamp2018-12-03', #make sure this is correct
                          row.uid, str(row.pre_cell_id), str(row.post_cell_id)), axis=1)
#df.to_csv('average_fit_vclamp_2018_12_06.csv')  #comment out after using so dont overwrite something