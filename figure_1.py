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

parser = argparse.ArgumentParser()
# specify each connection for example figure as exp ID, pre_cell, post_cell
parser.add_argument('--tlx3', nargs=3, type=int)
parser.add_argument('--sim1', nargs=3, type=int)
parser.add_argument('--rorb', nargs=3, type=int)
parser.add_argument('--l23pyr', nargs=3, type=int)

args = vars(parser.parse_args(sys.argv[1:]))
all_expts = ExperimentList(cache='expts_cache.pkl')

grid = PlotGrid()
grid.set_shape(4, 3)
grid.show()
grid_positions = {'l23pyr': 0, 'rorb': 1, 'sim1': 2, 'tlx3': 3}
grid[0,0].setTitle(title='First Pulse')
grid[0,1].setTitle(title='50 Hz Train')
grid[0,2].setTitle(title='Exponential Deconvolution')

for connection_type, expt_id in args.items():
    if expt_id is not None:
        expt_ind = expt_id[0]
        pre_cell = expt_id[1]
        post_cell = expt_id[2]
        expt = all_expts[expt_ind]
        row = grid_positions[connection_type]
        grid[row,0].addLegend()

        analyzer = DynamicsAnalyzer(expt, pre_cell, post_cell)
        if len(analyzer.pulse_responses) == 0:
            raise Exception("No suitable data found for cell %d -> cell %d in expt %s" % (pre_cell, post_cell, expt_ind))

        amp_group = analyzer.amp_group
        stop_freq = 50
        tau = 15e-3
        lp = 1000
        sweep_list = []
        ind = []
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
        if len(sweep_list) > 5:
            avg_first_pulse = trace_mean(sweep_list)
            avg_first_pulse.t0 = 0
            grid[row,0].setLabels(left=('Vm', 'V'))
            grid[row,0].setLabels(bottom=('t', 'sec'))
            grid[row,0].plot(avg_first_pulse.time_values, avg_first_pulse.data, name=connection_type)
        else:
            print ("%s not enough sweeps for first pulse" % connection_type)
        train_responses = analyzer.train_responses
        for i, stim_params in enumerate(train_responses.keys()):
            if stim_params[0] == 50:
                if len(train_responses[stim_params][0]) != 0:
                    ind_group = train_responses[stim_params][0]
                    for j in range(len(ind_group)):
                        ind.append(ind_group.responses[j])
        if len(ind) > 5:
            ind_avg = trace_mean(ind)
            ind_avg.t0 = 0
            base = float_mode(ind_avg.data[:int(10e-3 / ind_avg.dt)])
            ind_dec = bessel_filter(exp_deconvolve(ind_avg, tau), lp)
            ind_dec.t0 = 0
            grid[row,1].setLabels(left=('Vm','V'))
            grid[row,1].setLabels(bottom=('t', 'sec'))
            grid[row,2].setLabels(bottom=('t', 'sec'))
            grid[row,1].plot(ind_avg.time_values, ind_avg.data - base)
            grid[row,2].plot(ind_dec.time_values, ind_dec.data)
            grid[row,0].setYLink(grid[row, 1])
        else:
            print ("%s not enough sweeps for trains" % connection_type)
            grid[row,0].setYLink(grid[0,1])