import pyqtgraph as pg
import numpy as np
import csv
import sys
import argparse
from multipatch_analysis.experiment_list import cached_experiments
from manuscript_figures import pulse_qc, get_response, get_amplitude, response_filter
from synapse_comparison import load_cache, summary_plot_pulse
from neuroanalysis.data import TraceList
from neuroanalysis.ui.plot_grid import PlotGrid
from multipatch_analysis.connection_detection import fit_psp
from rep_connections import ee_connections, human_connections, no_include, all_connections, ie_connections, ii_connections, ei_connections
from multipatch_analysis.synaptic_dynamics import DynamicsAnalyzer
import matplotlib.pyplot as mplt
from scipy import stats
import time
import pandas as pd
import time as sys_time

parser = argparse.ArgumentParser(description='Enter organism and type of connection you"d like to analyze ex: mouse ee (all mouse excitatory-'
                'excitatory). Alternatively enter a cre-type connection ex: sim1-sim1')
parser.add_argument('--organism', dest='organism', help='Select mouse or human')
parser.add_argument('--connection', dest='connection', help='Specify connections to analyze')
args = vars(parser.parse_args(sys.argv[1:]))

all_expts = cached_experiments()

if args['organism'] == 'mouse':
    calcium = 'high'
    age = '40-60'
    sweep_threshold = 3
    threshold = 0.03e-3
    connection = args['connection']
    if connection == 'ee':
        connection_types = ee_connections.keys()
    elif connection == 'ii':
        connection_types = ii_connections.keys()
    elif connection == 'ei':
        connection_types = ei_connections.keys()
    elif connection == 'ie':
        connection_types == ie_connections.keys()
    elif connection == 'all':
        connection_types = all_connections.keys()
    elif len(connection.split('-')) == 2:
        conn_type = connection.split('-')
        if conn_type[0] == '2/3':
            pre_type = ('2/3', 'unknown')
        else:
            pre_type = (None, conn_type[0])
        if conn_type[1] == '2/3':
            post_type = ('2/3', 'unknown')
        else:
            post_type = (None, conn_type[0])
        connection_types = [(pre_type, post_type)]
elif args['organism'] == 'human':
    calcium = None
    age = None
    sweep_threshold = 5
    threshold = None
    connection = args['connection']
    if connection == 'ee':
        connection_types = human_connections.keys()
    else:
        conn_type = connection.split('-')
        connection_types = [((conn_type[0], 'unknown'), (conn_type[1], 'unknown'))]

holding = [-65, -75]
expt_ids = {}

for c in range(len(connection_types)):
    cre_type = (connection_types[c][0][1], connection_types[c][1][1])
#    if cre_type[0]!='rorb' and cre_type[1]!='rorb':
#        print 'skipping'
#        continue
    target_layer = (connection_types[c][0][0], connection_types[c][1][0])
    conn_type = connection_types[c]
#---here she does a selection I don't think I need to do    
    expt_list = all_expts.select(cre_type=cre_type, target_layer=target_layer, calcium=calcium, age=age)
#    grand_response[conn_type[0]] = {'trace': [], 'amp': [], 'latency': [], 'rise': [], 'dist': [], 'decay':[], 'CV': []}
    for expt in expt_list:
        for pre, post in expt.connections:
            if [expt.uid, pre, post] in no_include:
                continue
##            if str(expt.uid)!='1502301827.80':
#            if str(expt.uid)!='1501627688.56':
#            if str(expt.uid)!='1512537689.69': 
#                continue
            
            cre_check = expt.cells[pre].cre_type == cre_type[0] and expt.cells[post].cre_type == cre_type[1]
            layer_check = expt.cells[pre].target_layer == target_layer[0] and expt.cells[post].target_layer == target_layer[1]
            if cre_check is True and layer_check is True:
# ----Get response from pulse_response table
                # this gets response via the dynamics analyzer
                pulse_response, artifact = get_response(expt, pre, post, type='pulse')

                if threshold is not None and artifact > threshold:
                    continue
                response_subset, hold = response_filter(pulse_response, freq_range=[0, 50], holding_range=holding, pulse=True)
                if len(response_subset) >= sweep_threshold:
#                    qc_plot.clear()
                    qc_list = pulse_qc(response_subset, baseline=1.5, pulse=None, plot=None)
                    if len(qc_list) >= sweep_threshold:
#----here it is getting the amplitude so perhaps we don't need the response because the amplitude is cached                        
                        avg_trace, avg_amp, amp_sign, peak_t = get_amplitude(qc_list)
                        if amp_sign is '-':
                            continue
                        #print ('%s, %0.0f' %((expt.uid, pre, post), hold, ))
#                        all_amps = fail_rate(response_subset, '+', peak_t)
#                        cv = np.std(all_amps)/np.mean(all_amps)

                        psp_fits = fit_psp(avg_trace, sign=amp_sign, yoffset=0, amp=avg_amp, method='leastsq', stacked = False, fit_kws={})

##                        avg_trace.t0 = -(psp_fits.best_values['xoffset'] - 10e-3) #TODO: what is this
#                        expt_ids[conn_type[0]].append((pre, post, expt.uid, expt.source_id))
#

                        dt = avg_trace.dt
                        weight = np.ones(len(avg_trace.data))*10.
                        #weight[:int(10e-3/dt)] = 0.5
                        weight[int(10e-3/dt):int(12e-3/dt)] = 0  #area around stim artifact
                        weight[int(12e-3/dt):int(19e-3/dt)] = 30

                        time=avg_trace.time_values*1.e3

                        fig=mplt.figure(figsize=(20,8))
                        ax=fig.add_subplot(1,2,1)
                        ax2=ax.twinx()
                        ax.plot(time, avg_trace.data*1.e3)
                        ax2.plot(time, weight)
                        mplt.show()


##                        latency=(psp_fits.best_values['xoffset'] - 10e-3)*1.e3
##                        ax.plot(time, psp_fits.eval()*1.e3, label='NStack, factor=2, error='+str(np.around(psp_fits.nrmse(), 4))+', rise='+str(np.around(psp_fits.best_values['rise_time']*1e3, 4))+', latency='+str(np.around(latency, 4)))
##                        
##                        psp_fits = fit_psp(avg_trace, sign=amp_sign, yoffset=0, amp=avg_amp, method='leastsq', stacked = True, rise_time_mult_factor=2, fit_kws={})
##                        latency=(psp_fits.best_values['xoffset'] - 10e-3)*1.e3
##                        ax.plot(time, psp_fits.eval()*1.e3, label='Stack, factor=2, error='+str(np.around(psp_fits.nrmse(),4))+', rise='+str(np.around(psp_fits.best_values['rise_time']*1e3, 4))+', latency='+str(np.around(latency, 4)))
##                        
##                        psp_fits = fit_psp(avg_trace, sign=amp_sign, yoffset=0, amp=avg_amp, method='leastsq', stacked = True, rise_time_mult_factor=3, fit_kws={})
##                        latency=(psp_fits.best_values['xoffset'] - 10e-3)*1.e3
##                        ax.plot(time, psp_fits.eval()*1.e3, label='Stack, factor=3, error='+str(np.around(psp_fits.nrmse(), 4))+', rise='+str(np.around(psp_fits.best_values['rise_time']*1e3, 4))+', latency='+str(np.around(latency, 4)))                        
##                        
##                        psp_fits = fit_psp(avg_trace, sign=amp_sign, yoffset=0, amp=avg_amp, method='leastsq', stacked = False, rise_time_mult_factor=3, fit_kws={})
##                        latency=(psp_fits.best_values['xoffset'] - 10e-3)*1.e3
##                        ax.plot(time, psp_fits.eval()*1.e3, label='NStack, factor=3, error='+str(np.around(psp_fits.nrmse(),4))+', rise='+str(np.around(psp_fits.best_values['rise_time']*1e3, 4))+', latency='+str(np.around(latency, 4)))         
#
#                        print 'id', expt.uid, pre, post                       
# 
##                        psp_fits = fit_psp(avg_trace, sign=amp_sign, yoffset=0, xoffset=11e-3, amp=avg_amp, method='leastsq', stacked = False, rise_time_mult_factor=3, fit_kws={'weights': weight})
##                        latency=(psp_fits.best_values['xoffset'] - 10e-3)*1.e3
##                        ax.plot(time, psp_fits.eval()*1.e3, lw=8, label='NStack, xoff=11e-3, factor=3, weighted, error='+str(np.around(psp_fits.nrmse(),4))+', rise='+str(np.around(psp_fits.best_values['rise_time']*1e3, 4))+', latency='+str(np.around(latency, 4)))         
##
##                        if str(expt.uid)=='1494369390.68' and pre==3 and post==2:
##                            pass
##                        else:
##                            psp_fits = fit_psp(avg_trace, sign=amp_sign, yoffset=0, xoffset=11e-3, amp=avg_amp, method='leastsq', stacked = True, rise_time_mult_factor=3, fit_kws={'weights': weight})
##                            latency=(psp_fits.best_values['xoffset'] - 10e-3)*1.e3
##                            ax.plot(time, psp_fits.eval()*1.e3, lw=6, label='Stack, xoff=11e-3, factor=3, weighted, error='+str(np.around(psp_fits.nrmse(),4))+', rise='+str(np.around(psp_fits.best_values['rise_time']*1e3, 4))+', latency='+str(np.around(latency, 4)))         
##
##                        psp_fits = fit_psp(avg_trace, sign=amp_sign, yoffset=0, xoffset=14e-3, amp=avg_amp, method='leastsq', stacked = False, rise_time_mult_factor=3, fit_kws={'weights': weight})
##                        latency=(psp_fits.best_values['xoffset'] - 10e-3)*1.e3
##                        ax.plot(time, psp_fits.eval()*1.e3, lw=4, label='NStack, xoff=14e-3, factor=3, weighted, error='+str(np.around(psp_fits.nrmse(),4))+', rise='+str(np.around(psp_fits.best_values['rise_time']*1e3, 4))+', latency='+str(np.around(latency, 4)))         
##
##                        psp_fits = fit_psp(avg_trace, sign=amp_sign, yoffset=0, xoffset=14e-3, amp=avg_amp, method='leastsq', stacked = True, rise_time_mult_factor=3, fit_kws={'weights': weight})
##                        latency=(psp_fits.best_values['xoffset'] - 10e-3)*1.e3
##                        ax.plot(time, psp_fits.eval()*1.e3, lw=2, label='Stack, xoff=14e-3, factor=3, weighted, error='+str(np.around(psp_fits.nrmse(),4))+', rise='+str(np.around(psp_fits.best_values['rise_time']*1e3, 4))+', latency='+str(np.around(latency, 4)))         
#
#                        def try_init_conditions(init_xoffset):
#                            psp_fits={}
#                            for xoffset in init_xoffset:
#                                  
#                                try:
#                                    psp_fits[xoffset] = fit_psp(avg_trace, 
#                                                                sign=amp_sign, 
#                                                                yoffset=0, 
#                                                                xoffset=xoffset, 
#                                                                amp=avg_amp, 
#                                                                method='leastsq', 
#                                                                stacked = True, 
#                                                                rise_time_mult_factor=10., 
#                                                                fit_kws={'weights': weight})
#                                except:
#                                    psp_fits[xoffset] = 'fitting error'
#                        
#                            return psp_fits
#                        
#                        #fit average trace at a variety of initial conditions
#                        init_xoff_try=[11.e-3, 12.5e-3, 14.e-3]
#                        psp_fits=try_init_conditions(init_xoff_try)
#                        error=[]
#                        for xoffset in init_xoff_try:
#                            if xoffset not in psp_fits.keys():
#                                raise Exception('the offset should be in the keys')
#                            if psp_fits[xoffset]!='fitting error':
#                                error.append(psp_fits[xoffset].nrmse())
#                            else:
#                                error.append(np.infty)
#                        
#                        #get ind of min error
#                        min_index=np.where(np.array(error)==np.min(np.array(error)))[0][0] 
#                        best_key=init_xoff_try[min_index]
#                        for jj,key in enumerate(init_xoff_try):
#                            if psp_fits[key]!='fitting error':
#                                #make solid line the best fit
#                                latency=(psp_fits[key].best_values['xoffset'] - 10e-3)*1.e3
#                                nrmse=psp_fits[key].nrmse()
#                                rise_time=psp_fits[key].best_values['rise_time']*1e3
#                                if key==best_key:
#                                    ax.plot(time, 
#                                            psp_fits[key].eval()*1.e3, 
#                                            lw=9-jj*3, 
#                                            label='Stack, xoff='+str(key*1.e3)+
#                                                ', factor=10, weighted, error='+
#                                                str(np.around(nrmse,4))+
#                                                ', rise='+str(np.around(rise_time, 4))+
#                                                ', latency='+str(np.around(latency, 4)))         
#                                if key!=best_key:
#                                    ax.plot(time, 
#                                            psp_fits[key].eval()*1.e3, 
#                                            '--',
#                                            lw=9-jj*3, 
#                                            label='Stack, xoff='+str(key*1.e3)+
#                                                ', factor=10, weighted, error='+
#                                                str(np.around(nrmse,4))+
#                                                ', rise='+str(np.around(rise_time, 4))+
#                                                ', latency='+str(np.around(latency, 4)))                               
#
#                                                  
#                        mplt.title(str(len(qc_list))+' '+str(expt.uid)+' '+str(pre)+' '+str(post)+' '+str(conn_type))
#                        ax2.set_ylabel('weight')
#                        ax.set_ylabel('voltage (mV)')
#                        ax.set_xlabel('time (ms)')
#                        ax.legend(bbox_to_anchor=(1.05, 1), loc=2)
#                        ax.set_xlim([time[0], time[-1]])
#                        name=str(expt.uid)+'_'+str(pre)+'_'+str(post)
##                        mplt.show()
#                        mplt.savefig('/home/corinnet/Desktop/SynPhys_images/fitting/'+name+'.png')
#                        mplt.close()
#
#                    
##                    decay_response = response_filter(pulse_response, freq_range=[0, 20], holding_range=holding)
##                    qc_list = pulse_qc(response_subset, baseline=2, pulse=None)#, plot=qc_plot)
##                    if len(qc_list) >= sweep_threshold:
##                        avg_trace, avg_amp, amp_sign, peak_t = get_amplitude(qc_list)
##                        if amp_sign is '-':
##                            continue
##                        psp_fits = fit_psp(avg_trace, sign=amp_sign, yoffset=0, amp=avg_amp, method='leastsq', stacked = True,  fit_kws={})

##                        mplt.figure()
##                        mplt.plot(psp_fits.data)
##                        mplt.plot(psp_fits.eval())
##                        mplt.title(str(len(qc_list))+' '+str(expt.uid)+' '+str(pre)+' '+str(post)+' '+str(conn_type)+' '+str(psp_fits.nrmse()))
##                        mplt.show()

##write_cache(expt_ids, 'pulse_expt_ids.pkl')
##write_cache(features, 'pulse_features_human.pkl')
#
#df = pd.DataFrame(data=manifest)
#df = df[['Type', 'Connection', 'amp', 'latency', 'rise', 'decay', 'CV', 'nrmse']]
#writer = pd.ExcelWriter('Fig1_manifest.xlsx')
#df.to_excel(writer, 'Sheet1')
#writer.save()
