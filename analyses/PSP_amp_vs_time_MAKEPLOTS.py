import matplotlib.pyplot as plt
import allensdk.core.json_utilities as ju
import numpy as np


def convert_none_to_nan(data_list):
    '''convert list of data with 'None' in it to nan
    '''
    new_data_list=[np.nan if value==None else value for value in data_list]
    return new_data_list

dictionary=ju.read("PSP_vs_time_output_data/psp_vs_time_2mMCa.json")
for key0 in dictionary.keys():
    for key1 in dictionary[key0].keys():
        if key1!='num_of_synapses':
            dictionary[key0][key1]=convert_none_to_nan(dictionary[key0][key1])
    
for key in dictionary.keys():   
    time_in_min=[v/60. for v in dictionary[key]['time_points']]
    plt.errorbar(time_in_min, dictionary[key]['avg_data'],  yerr=dictionary[key]['std_err'], label=key+', n='+str(dictionary[key]['num_of_synapses']))

plt.title('Summary of deflection amplitute of first pulse over the\n course of an experiment (averaged across synapses)')
plt.legend(loc=4)
plt.ylabel('voltage (mV)')
plt.xlabel('time since first synaptic epoch in experiment (m)')
plt.show(block=False)
plt.show()