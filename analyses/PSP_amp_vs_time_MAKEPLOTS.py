import matplotlib.pyplot as plt
import allensdk.core.json_utilities as ju
import numpy as np
import scipy
import scipy.stats
import pickle


def convert_none_to_nan(data_list):
    '''convert list of data with 'None' in it to nan
    '''
    new_data_list=[np.nan if value==None else value for value in data_list]
    return new_data_list

#file=open('pulse_expt_ids.pkl')
#Stephs_data=pickle.load(file)
#Steph_uids=np.array([])
#for key in Stephs_data:
#    Steph_uids=np.append(Steph_uids, np.array([d[2] for d in Stephs_data[key]]))

dictionary=ju.read("PSP_vs_time_output_data/uid_included_02_20_18.json")
for syn_type in dictionary.keys():
    for key1 in dictionary[syn_type].keys():
        if key1 != 'num_of_synapses' and key1 != 'slopes' and key1 != 'intercepts':
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
for key3 in ['sweep_numbers', 'individual_start_times']:
    for syn_qc in ['correct_amp_good_HP', 'correct_amp_bad_HP', 'wrong_amp_good_HP','wrong_amp_bad_HP']:
        box_plot_data=[]
        plt.figure()
        for syn_type in np.sort(dictionary.keys()):
            if key3=='individual_start_times':
                x=np.array(dictionary[syn_type]['slopes'][key3][syn_qc])*1e6*60.  #converting to nV/min
                y_units=' (uV/min)'
            elif key3=='sweep_numbers':
                x=np.array(dictionary[syn_type]['slopes'][key3][syn_qc])*1e6  #converting to mV/sweep
                y_units=' (uV/sweep)'
            else:
                raise Exception('the y units are not a known type')
            box_plot_data.append(x)
        bp=plt.boxplot(box_plot_data, notch=0, whis=[5, 95], sym='.')
        for ii in range(len(box_plot_data)):
            plt.plot(np.ones(len(box_plot_data[ii]))+ii, box_plot_data[ii], '.k', ms=16)
        plt.xlim(0.5, 5.5)
        plt.xticks([1, 2,3,4, 5], np.sort(dictionary.keys()), rotation=40)
        plt.ylabel('slope of linear regression' + y_units)
        plt.title('Non-normalized to starting PSP amp\n'+key3+', '+syn_qc)
        plt.setp(bp['boxes'], color=colors[syn_qc], linewidth=3)
        plt.setp(bp['whiskers'], color=colors[syn_qc], linewidth=3)
        plt.setp(bp['fliers'], mfc=colors[syn_qc], markersize=12)
        plt.setp(bp['medians'], color=colors[syn_qc], linewidth=3)
        plt.setp(bp['caps'], color=colors[syn_qc], linewidth=3)
        plt.tight_layout()

#------------------------------------------------------------------
#---------compare regression PSP start to real PSP amp-------------
#------------------------------------------------------------------

# calculate PSP amp start via regression
for key3 in ['individual_start_times']:
    plt.figure()
    for syn_type in np.sort(dictionary.keys()):
        intercepts=np.array(dictionary[syn_type]['intercepts'][key3]['correct_amp_good_HP'])
        plt.plot(np.array(dictionary[syn_type]['PSPs_amp_start'])*1e3, intercepts*1e3, '.', MS=12, label=syn_type)
    plt.legend(loc=4)
    plt.xlabel('Experimental PSP amp at start (mV)')
    plt.ylabel('PSP amp at start via linear regression (mV)')
    
#------------------------------------------------------------------
#---compare regression slope versus starting PSP amplitude---------
#------------------------------------------------------------------

# calculate PSP amp start via regression
for key3 in ['individual_start_times']:
    plt.figure()
    for syn_type in np.sort(dictionary.keys()):
        intercepts=np.array(dictionary[syn_type]['intercepts'][key3]['correct_amp_good_HP'])
        slopes=np.array(dictionary[syn_type]['slopes'][key3]['correct_amp_good_HP'])
        plt.plot(slopes*1e6*60, intercepts*1e3, '.', MS=12, label=syn_type)
    plt.legend(loc=3)
    plt.xlabel('slopes (uV/min)')
    plt.ylabel('PSP amp at start via linear regression (mV)')
    
#------------------------------------------------------------------
#---------compare regression PSP at end to real PSP amp-------------
#------------------------------------------------------------------

for key3 in ['individual_start_times']:
    plt.figure()
    for syn_type in np.sort(dictionary.keys()):
        intercepts=np.array(dictionary[syn_type]['intercepts'][key3]['correct_amp_good_HP'])
        slopes=np.array(dictionary[syn_type]['slopes'][key3]['correct_amp_good_HP'])   
        time_end=np.array(dictionary[syn_type]['length_of_experiment'])
#        plt.figure()
        end_amps=slopes*time_end+intercepts # calculate PSP amp end via regression
#        plt.hist(time_end, label=syn_type, alpha=.5) #checked this it looks good
#        plt.hist(np.array(dictionary[syn_type]['PSPs_amp_end'])*1e3, label=syn_type, alpha=.5)
        plt.plot(np.array(dictionary[syn_type]['PSPs_amp_end'])*1e3, end_amps*1e3, '.', MS=12, label=syn_type)
        plt.legend(loc=4)
    plt.xlabel('Experimental PSP amp at end (mV)')
    plt.ylabel('PSP amp at end via linear regression (mV)')


#------------------------------------------------------------------
#---convert non normalized slopes to percent decrease-------------
#------------------------------------------------------------------

# divide slope by initial PSP height to get % psp decrease per min
for key3 in ['individual_start_times']:
    for syn_qc in ['correct_amp_good_HP']:
        box_plot_data=[]
        plt.figure()
        for syn_type in np.sort(dictionary.keys()):
            if key3=='individual_start_times':
                intercepts=np.array(dictionary[syn_type]['intercepts'][key3][syn_qc])
                slopes=np.array(dictionary[syn_type]['slopes'][key3][syn_qc])   
                time_end=np.array(dictionary[syn_type]['length_of_experiment'])
                end_amps=slopes*time_end+intercepts
                same_sign=np.divide(intercepts, end_amps)>0 #check if numbers are the same sign
                x=(slopes/intercepts)[same_sign]*60*100. #percent decrease per minute             
                y_units=' (%/min)'
            elif key3=='sweep_numbers':
                pass
            else:
                raise Exception('the y units are not a known type')
            box_plot_data.append(x)
        bp=plt.boxplot(box_plot_data, notch=0, whis=[5, 95], sym='.')
        for ii in range(len(box_plot_data)):
            plt.plot(np.ones(len(box_plot_data[ii]))+ii, box_plot_data[ii], '.k', ms=16)
        plt.xlim(0.5, 5.5)
        plt.xticks([1, 2,3,4, 5], np.sort(dictionary.keys()), rotation=40)
        plt.ylabel('PSP change over time' + y_units)
        plt.title('% decrease in PSP amplitude')
        plt.setp(bp['boxes'], color=colors[syn_qc], linewidth=3)
        plt.setp(bp['whiskers'], color=colors[syn_qc], linewidth=3)
        plt.setp(bp['fliers'], mfc=colors[syn_qc], markersize=12)
        plt.setp(bp['medians'], color=colors[syn_qc], linewidth=3)
        plt.setp(bp['caps'], color=colors[syn_qc], linewidth=3)
        plt.tight_layout()

#------------------------------------------------------------------
#---total percent decrease over the course of the experiment-------
#------------------------------------------------------------------

for key3 in ['individual_start_times']:
    for syn_qc in ['correct_amp_good_HP']:
        box_plot_data=[]
        plt.figure()
        for syn_type in np.sort(dictionary.keys()):
            if key3=='individual_start_times':
                intercepts=np.array(dictionary[syn_type]['intercepts'][key3][syn_qc])
                slopes=np.array(dictionary[syn_type]['slopes'][key3][syn_qc])   
                time_end=np.array(dictionary[syn_type]['length_of_experiment'])
                end_amps=slopes*time_end+intercepts
                same_sign=np.divide(intercepts, end_amps)>0 #check if numbers are the same sign
#                PSP_end_experimental=np.array(dictionary[syn_type]['PSPs_amp_end'])
#                PSP_start_experimental=np.array(dictionary[syn_type]['PSPs_amp_start'])
                x=np.divide((end_amps-intercepts),intercepts)[same_sign]          
                y_units=' (%)'
            elif key3=='sweep_numbers':
                pass
            else:
                raise Exception('the y units are not a known type')
            box_plot_data.append(x*100.)
            print syn_type, np.median(x*100.), np.mean(x*100), scipy.stats.sem(x*100)
        bp=plt.boxplot(box_plot_data, notch=0, whis=[5, 95], sym='.')
        for ii in range(len(box_plot_data)):
            plt.plot(np.ones(len(box_plot_data[ii]))+ii, box_plot_data[ii], '.k', ms=16)
        plt.xlim(0.5, 5.5)
        plt.xticks([1, 2,3,4, 5], np.sort(dictionary.keys()), rotation=40)
        plt.ylabel('Total loss in signal over time of experiment' + y_units)
        plt.setp(bp['boxes'], color=colors[syn_qc], linewidth=3)
        plt.setp(bp['whiskers'], color=colors[syn_qc], linewidth=3)
        plt.setp(bp['fliers'], mfc=colors[syn_qc], markersize=12)
        plt.setp(bp['medians'], color=colors[syn_qc], linewidth=3)
        plt.setp(bp['caps'], color=colors[syn_qc], linewidth=3)
        plt.tight_layout()
        print syn_type, np.median
        
plt.show()