import sys
import argparse
import pyqtgraph as pg
from experiment_list import ExperimentList
from synaptic_dynamics import DynamicsAnalyzer
from connection_detection import trace_mean
from neuroanalysis.baseline import float_mode
from neuroanalysis.filter import bessel_filter
from neuroanalysis.event_detection import exp_deconvolve
from neuroanalysis.ui.plot_grid import PlotGrid


app = pg.mkQApp()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
spike_color = (169, 169, 169)
stim_color = (211, 211, 211)
response_color = 'k'

parser = argparse.ArgumentParser()
# specify each connection for example figure as exp ID, pre_cell, post_cell
parser.add_argument('--tlx3', nargs=3, type=int)
parser.add_argument('--sim1', nargs=3, type=int)
parser.add_argument('--rorb', nargs=3, type=int)
parser.add_argument('--l23pyr', nargs=3, type=int)

args = vars(parser.parse_args(sys.argv[1:]))
all_expts = ExperimentList(cache='expts_cache.pkl')

grid = PlotGrid()
grid.set_shape(9, 3)
grid.show()
for g in range(1, grid.shape[0], 2):
    grid.grid.ci.layout.setRowStretchFactor(g, 5)
    grid.grid.ci.layout.setRowStretchFactor(g + 1, 20)
    grid[g, 0].hideAxis('bottom')
    grid[g, 1].hideAxis('bottom')
    grid[g, 2].hideAxis('bottom')
    grid[g, 2].hideAxis('left')
grid.grid.ci.layout.setRowStretchFactor(0, 1)
grid[0, 0].hideAxis('bottom')
grid[0, 1].hideAxis('bottom')
grid[0, 2].hideAxis('bottom')
grid[0, 2].hideAxis('left')

grid_positions = {'l23pyr': range(1, 3), 'rorb': range(3, 5), 'sim1': range(5, 7), 'tlx3': range(7, 9)}
grid[0,0].setTitle(title='First Pulse')
grid[0,1].setTitle(title='50 Hz Train')
grid[0,2].setTitle(title='Exponential Deconvolution')
ii = 0
for connection_type, expt_id in args.items():
    if expt_id is not None:
        ii += 1
        expt_ind = expt_id[0]
        pre_cell = expt_id[1]
        post_cell = expt_id[2]
        expt = all_expts[expt_ind]
        row = grid_positions[connection_type]
        grid[row[1],0].addLegend()

        analyzer = DynamicsAnalyzer(expt, pre_cell, post_cell)
        if len(analyzer.pulse_responses) == 0:
            raise Exception("No suitable data found for cell %d -> cell %d in expt %s" % (pre_cell, post_cell, expt_ind))

        amp_group = analyzer.amp_group
        stop_freq = 50
        tau = 15e-3
        lp = 1000
        sweep_list = {'response': [], 'spike': [], 'command': []}
        ind = {'response': [], 'spike': []}
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
                post_base = float_mode(sweep_trace.data[:int(10e-3 / sweep_trace.dt)])
                pre_spike = amp_group.spikes[sweep]
                pre_base = float_mode(pre_spike.data[:int(10e-3 / pre_spike.dt)])
                sweep_list['response'].append(sweep_trace.copy(data=sweep_trace.data - post_base))
                sweep_list['spike'].append(pre_spike.copy(data=pre_spike.data - pre_base))
                sweep_list['command'].append(amp_group.commands[sweep])
        if len(sweep_list['response']) > 5:
            avg_first_pulse = trace_mean(sweep_list['response'])
            avg_first_pulse.t0 = 0
            avg_spike = trace_mean(sweep_list['spike'])
            avg_spike.t0 = 0
            grid[row[1],0].setLabels(left=('Vm', 'V'))
            grid[row[1],0].setLabels(bottom=('t', 's'))
            grid[row[1],0].plot(avg_first_pulse.time_values, avg_first_pulse.data, pen={'color': response_color, 'width': 2}, name=connection_type)
            grid[row[0],0].setLabels(left=('Vm', 'V'))
            grid[row[0],0].plot(avg_spike.time_values, avg_spike.data, pen=spike_color)
        else:
            print ("%s not enough sweeps for first pulse" % connection_type)
        if ii == 1:
            stim_command = trace_mean(sweep_list['command'])
            stim_command.t0 = 0
            grid[0, 0].plot(stim_command.time_values, stim_command.data, pen=stim_color)
            grid[0, 0].setLabels(left=('', ''))
            grid[0, 0].getAxis('left').setOpacity(0)
        train_responses = analyzer.train_responses
        for i, stim_params in enumerate(train_responses.keys()):
            if stim_params[0] == 50:
                if len(train_responses[stim_params][0]) != 0:
                    ind_group = train_responses[stim_params][0]
                    for j in range(len(ind_group)):
                        ind['response'].append(ind_group.responses[j])
                        ind['spike'].append(ind_group.spikes[j])
        if len(ind['response']) > 5:
            ind_avg = trace_mean(ind['response'])
            ind_avg.t0 = 0
            ind_base = float_mode(ind_avg.data[:int(10e-3 / ind_avg.dt)])
            ind_dec = bessel_filter(exp_deconvolve(ind_avg, tau), lp)
            ind_dec.t0 = 0
            train_spike = trace_mean(ind['spike'])
            train_spike.t0 = 0
            spike_base = float_mode(train_spike.data[:int(10e-3 / train_spike.dt)])
            stim_command = trace_mean(ind_group.commands)
            stim_command.t0 = 0
            grid[row[1],1].setLabels(left=('Vm','V'))
            grid[row[1],1].setLabels(bottom=('t', 's'))
            grid[row[1],2].setLabels(bottom=('t', 's'))
            grid[row[1],1].plot(ind_avg.time_values, ind_avg.data - ind_base, pen={'color': response_color, 'width': 2})
            grid[row[0], 1].setLabels(left=('Vm', 'V'))
            grid[row[0], 1].plot(train_spike.time_values, train_spike.data - spike_base, pen=spike_color)
            grid[row[1],2].plot(ind_dec.time_values, ind_dec.data, pen={'color': response_color, 'width': 2})
            grid[row[1],0].setYLink(grid[row[1], 1])
        else:
            print ("%s not enough sweeps for trains" % connection_type)
            grid[row[1], 0].setYLink(grid[2, 1])
            grid[row[0], 1].hideAxis('bottom')
            grid[row[0], 2].hideAxis('bottom')
            grid[row[0], 1].hideAxis('left')
            grid[row[0], 2].hideAxis('left')
            grid[row[1], 1].hideAxis('bottom')
            grid[row[1], 2].hideAxis('bottom')
            grid[row[1], 1].hideAxis('left')
            grid[row[1], 2].hideAxis('left')
        if ii == 1:
            stim_command = trace_mean(ind_group.commands)
            stim_command.t0 = 0
            grid[0, 1].plot(stim_command.time_values, stim_command.data, pen=stim_color)
            grid[0, 1].setLabels(left=('', ''))
            grid[0, 1].getAxis('left').setOpacity(0)