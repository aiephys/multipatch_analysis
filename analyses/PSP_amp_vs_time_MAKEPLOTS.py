import matplotlib.pyplot as plt
import allensdk.core.json_utilities as ju
import numpy as np


def convert_none_to_nan(data_list):
    '''convert list of data with 'None' in it to nan
    '''
    new_data_list=[np.nan if value==None else value for value in data_list]
    return new_data_list

dictionary=ju.read("PSP_vs_time_output_data/goodpsp_vs_time or_sweep_2_05_18.json")
for syn_type in dictionary.keys():
    for key1 in dictionary[syn_type].keys():
        if key1 != 'num_of_synapses' and key1 != 'slopes':
            print key1
            dictionary[syn_type][key1]=convert_none_to_nan(dictionary[syn_type][key1])

plt.figure()    
for syn_type in dictionary.keys():   
    time_in_min=[v/60. for v in dictionary[syn_type]['time_points']]
    plt.errorbar(time_in_min, dictionary[syn_type]['time_avg_data'],  
                 yerr=dictionary[syn_type]['time_std_err'], label=syn_type+', n='+str(dictionary[syn_type]['num_of_synapses']))

plt.title('Summary of deflection amplitude of first pulse over the\n course of an experiment (averaged across synapses)')
plt.legend(loc=1)
plt.ylabel('voltage (V)')
plt.xlabel('time since first synaptic epoch in experiment (m)')
plt.show(block=False)


plt.figure()    
for syn_type in dictionary.keys():   
    plt.errorbar(dictionary[syn_type]['sweeps'], dictionary[syn_type]['sweep_avg_data'],  
                 yerr=dictionary[syn_type]['sweep_std_err'], label=syn_type+', n='+str(dictionary[syn_type]['num_of_synapses']))

plt.title('Summary of deflection amplitude of first pulse over the\n course of an experiment (averaged across synapses)')
plt.legend(loc=4)
plt.ylabel('voltage (V)')
plt.xlabel('sweep number')
plt.show(block=False)


#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------plotting distributions of slopes of the PSP amp over time or sweep number-----------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

def make_box_plot(data, the_color, p=[1,2,3,4,5]):
    bp=plt.boxplot([(data, 2),
                    (data, 3),
                    (data, 4),
                    (data, 5),
                    (data, 6)], notch=0, whis=[5, 95], sym='.', positions=p, widths=.08)
    plt.setp(bp['boxes'], color=the_color, linewidth=3)
    plt.setp(bp['whiskers'], color=the_color, linewidth=3)
    plt.setp(bp['fliers'], mfc=the_color, markersize=12)
    plt.setp(bp['medians'], color=the_color, linewidth=3)
    plt.setp(bp['caps'], color=the_color, linewidth=3)


colors={}
colors['correct_amp_good_HP']='g'
colors['correct_amp_bad_HP']='purple'
colors['wrong_amp_good_HP']='c'
colors['wrong_amp_bad_HP']='r'
for key3 in ['sweep_number', 'time']:
    for syn_qc in ['correct_amp_good_HP', 'correct_amp_bad_HP', 'wrong_amp_good_HP','wrong_amp_bad_HP']:
        box_plot_data=[]
        plt.figure()
        for syn_type in np.sort(dictionary.keys()):
            if key3=='time':
                x=np.array(dictionary[syn_type]['slopes'][key3][syn_qc])*1e6*60.  #converting to nV/min
                y_units=' (uV/min)'
            elif key3=='sweep_number':
                x=np.array(dictionary[syn_type]['slopes'][key3][syn_qc])*1e6  #converting to mV/sweep
                y_units=' (uV/sweep)'
            else:
                raise Exception('the y units are not a known type')
            box_plot_data.append(x)
        bp=plt.boxplot(box_plot_data, notch=0, whis=[5, 95], sym='.')
        plt.xlim(0.5, 5.5)
        plt.xticks([1, 2,3,4, 5], np.sort(dictionary.keys()), rotation=40)
        plt.ylabel('slope of linear regression' + y_units)
        plt.title(key3+', '+syn_qc)
        plt.setp(bp['boxes'], color=colors[syn_qc], linewidth=3)
        plt.setp(bp['whiskers'], color=colors[syn_qc], linewidth=3)
        plt.setp(bp['fliers'], mfc=colors[syn_qc], markersize=12)
        plt.setp(bp['medians'], color=colors[syn_qc], linewidth=3)
        plt.setp(bp['caps'], color=colors[syn_qc], linewidth=3)
        plt.tight_layout()

    

plt.show()

