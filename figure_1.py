import sys
import argparse
import pyqtgraph as pg
from experiment_list import ExperimentList
from synaptic_dynamics import DynamicsAnalyzer
from neuroanalysis.data import TraceList
from neuroanalysis.baseline import float_mode
from neuroanalysis.filter import bessel_filter
from neuroanalysis.event_detection import exp_deconvolve
from neuroanalysis.ui.plot_grid import PlotGrid


app = pg.mkQApp()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
grey = (169, 169, 169)

parser = argparse.ArgumentParser()
# specify each connection for example figure as exp ID, pre_cell, post_cell
parser.add_argument('--tlx3', nargs=3, type=int)
parser.add_argument('--sim1', nargs=3, type=int)
parser.add_argument('--rorb', nargs=3, type=int)
parser.add_argument('--l23pyr', nargs=3, type=int)
parser.add_argument('--sweeps', action='store_true', default=False, dest='sweeps',
                    help='plot individual sweeps behing average')
parser.add_argument('--trains', action='store_true', default=False, dest='trains', help='plot 50Hz train and deconvolution')
parser.add_argument('--axis', type=str, help='choose linked (to link all grid y-axes) or overlay (generate new'
                                                      'plot with first pulse responses overlaid)')

args = vars(parser.parse_args(sys.argv[1:]))
plot_sweeps =args['sweeps']
plot_trains = args['trains']
axis = args['axis']
all_expts = ExperimentList(cache='expts_cache.pkl')

grid = PlotGrid()
if plot_trains is True:
    grid.set_shape(9, 3)
else:
    grid.set_shape(9,1)
grid.show()
for g in range(1, grid.shape[0], 2):
    grid.grid.ci.layout.setRowStretchFactor(g, 5)
    grid.grid.ci.layout.setRowStretchFactor(g + 1, 20)
    grid[g, 0].hideAxis('bottom')
    if plot_trains is True:
        grid[g, 1].hideAxis('bottom')
        grid[g, 2].hideAxis('bottom')
        grid[g, 2].hideAxis('left')
        grid[g + 1, 2].hideAxis('left')
grid.grid.ci.layout.setRowStretchFactor(0, 1)
grid[0, 0].hideAxis('bottom')
if plot_trains is True:
    grid[0, 1].hideAxis('bottom')
    grid[0, 2].hideAxis('bottom')
    grid[0, 2].hideAxis('left')
if axis == 'overlay':
    overlay_plot = pg.plot()
    overlay_plot.addLegend()

grid_positions = {'l23pyr': range(1, 3), 'rorb': range(3, 5), 'sim1': range(5, 7), 'tlx3': range(7, 9)}
grid[0, 0].setTitle(title='First Pulse')
if plot_trains is True:
    grid[0, 1].setTitle(title='50 Hz Train')
    grid[0, 2].setTitle(title='Exponential Deconvolution')
ii = 0
for connection_type in grid_positions.keys():
    if args[connection_type] is not None:
        ii += 1
        expt_ind = args[connection_type][0]
        pre_cell = args[connection_type][1]
        post_cell = args[connection_type][2]
        expt = all_expts[expt_ind]
        row = grid_positions[connection_type]

        analyzer = DynamicsAnalyzer(expt, pre_cell, post_cell)
        if len(analyzer.pulse_responses) == 0:
            raise Exception("No suitable data found for cell %d -> cell %d in expt %s" % (pre_cell, post_cell, expt_ind))

        amp_group = analyzer.amp_group
        stop_freq = 50
        tau = 15e-3
        lp = 1000
        sweep_list = {'response': [], 'spike': [], 'command': []}
        ind = {'response': [], 'spike': [], 'dec': []}
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
            n = len(sweep_list['response'])
            if plot_sweeps is True:
                for sweep in range(n):
                    current_sweep = sweep_list['response'][sweep]
                    current_spike = sweep_list['spike'][sweep]
                    grid[row[1], 0].plot(current_sweep.time_values, current_sweep.data, pen=grey)
                    grid[row[0], 0].plot(current_spike.time_values, current_spike.data, pen=grey)
            avg_first_pulse = TraceList(sweep_list['response']).mean()
            avg_first_pulse.t0 = 0
            avg_spike = TraceList(sweep_list['spike']).mean()
            avg_spike.t0 = 0
            grid[row[1], 0].setLabels(left=('Vm', 'V'))
            grid[row[1], 0].setLabels(bottom=('t', 's'))
            grid[row[1], 0].setXRange(0, 30e-3)
            grid[row[1], 0].plot(avg_first_pulse.time_values, avg_first_pulse.data, pen={'color': 'k', 'width': 2})
            grid[row[0], 0].setLabels(left=('Vm', 'V'))
            grid[row[0], 0].plot(avg_spike.time_values, avg_spike.data, pen='k')
            grid[row[0], 0].setXLink(grid[row[1], 0])
            label = pg.LabelItem('%s, n = %d' % (connection_type, n))
            label.setParentItem(grid[row[1], 0].vb)
            label.setPos(50, 0)
            grid[row[1], 0].label = label
        else:
            print ("%s not enough sweeps for first pulse" % connection_type)
        if ii == 1:
            stim_command = TraceList(sweep_list['command']).mean()
            stim_command.t0 = 0
            grid[0, 0].plot(stim_command.time_values, stim_command.data, pen=grey)
            grid[0, 0].setLabels(left=('', ''))
            grid[0, 0].getAxis('left').setOpacity(0)
            grid[0, 0].setXLink(grid[2, 0])
        if axis == 'overlay':
            overlay_plot.plot(avg_first_pulse.time_values, avg_first_pulse.data, pen=(ii, len(grid_positions.keys())*1.3), name=connection_type)
            overlay_plot.setLabels(left=('Vm', 'V'))
            overlay_plot.setLabels(bottom=('t', 's'))
        if plot_trains is True:
            train_responses = analyzer.train_responses
            for i, stim_params in enumerate(train_responses.keys()):
                if stim_params[0] == 50:
                    if len(train_responses[stim_params][0]) != 0:
                        ind_group = train_responses[stim_params][0]
                        for j in range(len(ind_group)):
                            train_trace = ind_group.responses[j]
                            ind_base = float_mode(train_trace.data[:int(10e-3 / train_trace.dt)])
                            ind['response'].append(train_trace.copy(data=train_trace.data - ind_base))
                            dec_trace = bessel_filter(exp_deconvolve(train_trace, tau), lp)
                            dec_base = float_mode(dec_trace.data[:int(10e-3 / dec_trace.dt)])
                            ind['dec'].append(dec_trace.copy(data=dec_trace.data - dec_base))
                            ind['spike'].append(ind_group.spikes[j])
            if len(ind['response']) > 5:
                n = len(ind['response'])
                if plot_sweeps is True:
                    for sweep in range(n):
                        train_sweep = ind['response'][sweep]
                        dec_sweep = ind['dec'][sweep]
                        grid[row[1], 1].plot(train_sweep.time_values, train_sweep.data, pen=grey)
                        grid[row[1], 2].plot(dec_sweep.time_values, dec_sweep.data, pen=grey)
                ind_avg = TraceList(ind['response']).mean()
                ind_avg.t0 = 0
                ind_dec = TraceList(ind['dec']).mean()
                ind_dec.t0 = 0
                train_spike = TraceList(ind['spike']).mean()
                train_spike.t0 = 0
                spike_base = float_mode(train_spike.data[:int(10e-3 / train_spike.dt)])
                stim_command = TraceList(ind_group.commands).mean()
                stim_command.t0 = 0
                grid[row[1], 1].setLabels(left=('Vm','V'))
                grid[row[1], 1].setLabels(bottom=('t', 's'))
                grid[row[1], 2].setLabels(bottom=('t', 's'))
                grid[row[1], 1].plot(ind_avg.time_values, ind_avg.data, pen={'color': 'k', 'width': 2})
                grid[row[0], 1].setLabels(left=('Vm', 'V'))
                grid[row[0], 1].plot(train_spike.time_values, train_spike.data - spike_base, pen=grey)
                grid[row[1], 2].plot(ind_dec.time_values, ind_dec.data, pen={'color': 'k', 'width': 2})
                grid[row[1], 0].setYLink(grid[row[1], 1])
                label = pg.LabelItem('n = %d' % n)
                label.setParentItem(grid[row[1], 1].vb)
                label.setPos(50, 0)
                grid[row[1], 1].label = label
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
                stim_command = TraceList(ind_group.commands).mean()
                stim_command.t0 = 0
                grid[0, 1].plot(stim_command.time_values, stim_command.data, pen=grey)
                grid[0, 1].setLabels(left=('', ''))
                grid[0, 1].getAxis('left').setOpacity(0)
                grid[0, 1].setXLink(grid[2, 1])