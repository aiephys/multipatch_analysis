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
from synapse_comparison import load_cache
from manuscript_figures import cache_response, response_filter, trace_avg, pulse_qc, write_cache
from rep_connections import connections
from constants import INHIBITORY_CRE_TYPES, EXCITATORY_CRE_TYPES

### Select synapses for representative traces as {Connection Type: [UID, Pre_cell, Post_cell], } ###

connection_types = connections.keys() #['L23pyr-L23pyr', 'rorb-rorb', 'sim1-sim1', 'tlx3-tlx3']

cache_file = 'pulse_response_cache.pkl'
response_cache = load_cache(cache_file)
cache_change = []
pg.dbg()
app = pg.mkQApp()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
grey = (169, 169, 169)
sweep_color = (0, 0, 0, 30)
holding_i = [-53, -60]
holding_e = [-68, -72]

parser = argparse.ArgumentParser()
# plot options
parser.add_argument('--sweeps', action='store_true', default=False, dest='sweeps',
                    help='plot individual sweeps behing average')
parser.add_argument('--trains', action='store_true', default=False, dest='trains', help='plot 50Hz train and deconvolution')
parser.add_argument('--link-y-axis', action='store_true', default=False, dest='link-y-axis', help='link all y-axis down a column')

args = vars(parser.parse_args(sys.argv[1:]))
plot_sweeps = args['sweeps']
plot_trains = args['trains']
link_y_axis = args['link-y-axis']
all_expts = ExperimentList(cache='expts_cache.pkl')

grid = PlotGrid()
if plot_trains is True:
    grid.set_shape(len(connection_types), 3)
    train_cache_file = 'train_response_cache.pkl'
    train_response_cache = load_cache(train_cache_file)
    grid[0, 2].hideAxis('left')
    grid[0, 1].setTitle(title='50 Hz Train')
    grid[0, 2].setTitle(title='Exponential Deconvolution')
    tau = 15e-3
    lp = 1000
else:
    grid.set_shape(len(connection_types), 2)
grid.show()
row = 0
grid[0, 0].setTitle(title='First Pulse')

maxYpulse = []
maxYtrain = []
maxYdec = []
for row in range(len(connection_types)):
    connection_type = connection_types[row]
    if type(connection_type) is not tuple:
        connection_type = tuple(connection_type.split('-'))
    expt_id, pre_cell, post_cell = connections[connection_type]
    expt = all_expts[expt_id]
    if connection_type[0] in EXCITATORY_CRE_TYPES:
        holding = holding_e
        sign = '+'
    elif connection_type[0] in INHIBITORY_CRE_TYPES:
        holding = holding_i
        sign = '-'
    pulse_response, change = cache_response(expt, pre_cell, post_cell, response_cache, type='pulse')
    cache_change.append(change)
    sweep_list = response_filter(pulse_response, freq_range=[0, 50], holding_range=holding, pulse=True)
    n_sweeps = len(sweep_list)
    if n_sweeps > 5:
        qc_list = pulse_qc(sweep_list, baseline=4, pulse=4, plot=grid[row, 1])
        avg_first_pulse = trace_avg(sweep_list)
        avg_first_pulse.t0 = 0
        if plot_sweeps is True:
            for current_sweep in sweep_list:
                current_sweep.t0 = 0
                base = float_mode(current_sweep.data[:int(10e-3 / current_sweep.dt)])
                grid[row, 0].plot(current_sweep.time_values, current_sweep.data - base, pen=sweep_color)

        grid[row, 0].setLabels(left=('Vm', 'V'))
        grid[row, 0].setLabels(bottom=('t', 's'))
        grid[row, 0].setXRange(-2e-3, 27e-3)
        grid[row, 0].plot(avg_first_pulse.time_values, avg_first_pulse.data, pen={'color': (255, 0, 255), 'width': 2})
        grid[row, 0].setLabels(left=('Vm', 'V'))
        label = pg.LabelItem('%s, n = %d' % (connection_type, n_sweeps))
        label.setParentItem(grid[row, 0].vb)
        label.setPos(50, 0)
        maxYpulse.append((row, grid[row,0].getAxis('left').range[1]))
    else:
        print ("%s -> %s not enough sweeps for first pulse" % (connection_type[0], connection_type[1]))
    if plot_trains is True:
        train_responses = cache_response(expt, pre_cell, post_cell, train_cache_file, train_response_cache, type='train')
        train_sweep_list = response_filter(train_responses[0], freq_range=[50, 50], holding_range=holding)
        n_train_sweeps = len(train_sweep_list)
        if n_train_sweeps > 5:
            dec_sweep_list = []
            for sweep in range(n_train_sweeps):
                train_sweep = train_sweep_list[sweep]
                train_base = float_mode(train_sweep.data[:int(10e-3 / train_sweep.dt)])
                dec_sweep = bessel_filter(exp_deconvolve(train_sweep, tau), lp)
                dec_base = float_mode(dec_sweep.data[:int(10e-3 / dec_sweep.dt)])
                dec_sweep_list.append(dec_sweep.copy(data=dec_sweep.data - dec_base))
                if plot_sweeps is True:
                    grid[row, 1].plot(train_sweep.time_values, train_sweep.data - train_base, pen=sweep_color)
                    grid[row, 2].plot(dec_sweep.time_values, dec_sweep.data - dec_base, pen=sweep_color)
            ind_avg = trace_avg(train_sweep_list)
            ind_avg.t0 = 0
            ind_dec = trace_avg(dec_sweep_list)
            ind_dec.t0 = 0
            grid[row, 1].setLabels(left=('Vm','V'))
            grid[row, 1].setLabels(bottom=('t', 's'))
            grid[row, 2].setLabels(bottom=('t', 's'))
            grid[row, 1].plot(ind_avg.time_values, ind_avg.data, pen={'color': (255, 0, 255), 'width': 2})
            grid[row, 1].setLabels(left=('Vm', 'V'))
            grid[row, 2].plot(ind_dec.time_values, ind_dec.data, pen={'color': (255, 0, 255), 'width': 2})
            label = pg.LabelItem('n = %d' % n_train_sweeps)
            label.setParentItem(grid[row, 1].vb)
            label.setPos(50, 0)
            grid[row, 1].label = label
            maxYtrain.append((row, grid[row, 1].getAxis('left').range[1]))
            maxYdec.append((row, grid[row, 2].getAxis('left').range[1]))
        else:
            print ("%s not enough sweeps for trains" % connection_type)
            grid[row, 0].setYLink(grid[2, 1])
            grid[row, 1].hideAxis('bottom')
            grid[row, 2].hideAxis('bottom')
            grid[row, 1].hideAxis('left')
            grid[row, 2].hideAxis('left')
    row += 1
if link_y_axis is True:
    max_row = max(maxYpulse, key=lambda x:x[1])[0]
    for g in range(1, grid.shape[0], 2):
        if g != max_row:
            grid[g, 0].setYLink(grid[max_row, 0])
        if plot_trains is True:
            max_row_train = max(maxYtrain, key=lambda x:x[1])[0]
            max_row_dec = max(maxYdec, key=lambda x:x[1])[0]
            for g in range(2, grid.shape[0], 2):
                if g != max_row_train:
                    grid[g, 1].setYLink(grid[max_row_train, 1])
                if g != max_row_dec:
                    grid[g, 2].setYLink(grid[max_row_dec, 2])

# if sum(cache_change) > 0:
#     write_cache(response_cache, cache_file)

if sys.flags.interactive == 0:
    app.exec_()