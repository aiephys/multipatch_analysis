import matplotlib.pyplot as plt
import allensdk.core.json_utilities as ju
import numpy as np


def convert_none_to_nan(data_list):
    '''convert list of data with 'None' in it to nan
    '''
    new_data_list=[np.nan if value==None else value for value in data_list]
    return new_data_list

dictionary=ju.read("PSP_vs_time_output_data/goodpsp_vs_time or_sweep_1_29_18.json")
for key0 in dictionary.keys():
    for key1 in dictionary[key0].keys():
        if key1!='num_of_synapses':
            dictionary[key0][key1]=convert_none_to_nan(dictionary[key0][key1])

plt.figure()    
for key in dictionary.keys():   
    time_in_min=[v/60. for v in dictionary[key]['time_points']]
    plt.errorbar(time_in_min, dictionary[key]['time_avg_data'],  yerr=dictionary[key]['time_std_err'], label=key+', n='+str(dictionary[key]['num_of_synapses']))

plt.title('Summary of deflection amplitude of first pulse over the\n course of an experiment (averaged across synapses)')
plt.legend(loc=1)
plt.ylabel('voltage (V)')
plt.xlabel('time since first synaptic epoch in experiment (m)')
plt.show(block=False)

plt.figure()    
for key in dictionary.keys():   
    plt.errorbar(dictionary[key]['sweeps'], dictionary[key]['sweep_avg_data'],  yerr=dictionary[key]['sweep_std_err'], label=key+', n='+str(dictionary[key]['num_of_synapses']))

plt.title('Summary of deflection amplitude of first pulse over the\n course of an experiment (averaged across synapses)')
plt.legend(loc=4)
plt.ylabel('voltage (V)')
plt.xlabel('sweep number')
plt.show(block=False)
plt.show()