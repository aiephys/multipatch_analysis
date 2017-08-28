import sys
import pyqtgraph as pg
from experiment_list import ExperimentList
from synaptic_dynamics import DynamicsAnalyzer
from connection_detection import trace_mean
from neuroanalysis.baseline import float_mode


app = pg.mkQApp()
pg.dbg()

arg = sys.argv[1]
expt_ind = int(arg)
all_expts = ExperimentList(cache='expts_cache.pkl')
expt = all_expts[expt_ind]

pre_cell = int(sys.argv[2])
post_cell = int(sys.argv[3])

analyzer = DynamicsAnalyzer(expt, pre_cell, post_cell)
if len(analyzer.pulse_responses) == 0:
    raise Exception("No suitable data found for cell %d -> cell %d in expt %s" % (pre_cell, post_cell, expt_ind))

amp_group = analyzer.amp_group
stop_freq = 50
sweep_list = []
ind = []
ind_base_subtract = []
n_sweeps = len(amp_group)
if n_sweeps == 0:
    print "No Sweeps"
for sweep in range(n_sweeps):
    stim_name = amp_group.responses[sweep].recording.meta['stim_name']
    stim_param = stim_name.split('_')
    freq = stim_param[1]
    freq = int(freq.split('H')[0])
    if freq <= stop_freq:
        sweep_trace = amp_group.responses[sweep]
        base = float_mode(sweep_trace.data[:int(10e-3 / sweep_trace.dt)])
        sweep_list.append(sweep_trace.copy(data=sweep_trace.data - base))
if len(sweep_list) != 0:
    avg_first_pulse = trace_mean(sweep_list)
    amp_plot = pg.plot(title='First pulse amplitude')
    amp_plot.setLabels(left=('Vm', 'V'))
    amp_plot.plot(avg_first_pulse.time_values, avg_first_pulse.data)
train_responses = analyzer.train_responses
for i, stim_params in enumerate(train_responses.keys()):
    if stim_params[0] == 50:
        if len(train_responses[stim_params][0]) != 0:
            ind_group = train_responses[stim_params][0]
            for j in range(len(ind_group)):
                ind.append(ind_group.responses[j])
if len(ind) > 5:
    ind_avg = trace_mean(ind)
    base = float_mode(ind_avg.data[:int(10e-3 / ind_avg.dt)])
    ind_base_subtract.append(ind_avg.copy(data=ind_avg.data - base))
    train_plot = pg.plot(title='50 Hz train')
    train_plot.setLabels(left=('Vm','V'))
    train_plot.plot(ind_avg.time_values, ind_avg.data - base)




