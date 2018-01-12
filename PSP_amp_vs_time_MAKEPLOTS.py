import matplotlib.pyplot as plt
import allensdk.core.json_utilities as ju
import numpy as np


def convert_none_to_nan(data_list):
    '''convert list of data with 'None' in it to nan
    '''
    new_data_list=[np.nan if value==None else value for value in data_list]
    return new_data_list

dictionary=ju.read("PSP_vs_time_output_data/psp_vs_time1_12_18_R_drive_data.json")
for key0 in dictionary.keys():
    for key1 in dictionary[key0].keys():        
        dictionary[key0][key1]=convert_none_to_nan(dictionary[key0][key1])
    
for key in dictionary.keys():   
    plt.errorbar(dictionary[key]['time_points'], dictionary[key]['avg_data'],  yerr=dictionary[key]['std_err'], label=key)

plt.title('average base-line subtracted first pulse synaptic deflection')
plt.legend()
plt.ylabel('voltage (mV)')
plt.xlabel('time since first recorded synapse (s)')
plt.show(block=False)
plt.show()