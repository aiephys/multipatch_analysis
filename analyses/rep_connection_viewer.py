import sys
import argparse
import pyqtgraph as pg
import numpy as np
from aisynphys.experiment_list import cached_experiments
from neuroanalysis.filter import bessel_filter
from neuroanalysis.event_detection import exp_deconvolve
from neuroanalysis.ui.plot_grid import PlotGrid
from manuscript_figures import response_filter, trace_avg, pulse_qc, get_response, bsub, trace_plot, fail_rate, get_amplitude
from rep_connections import all_connections, human_connections, ee_connections
from aisynphys.constants import INHIBITORY_CRE_TYPES, EXCITATORY_CRE_TYPES

### Select synapses for representative traces as {('pre-type'), ('post-type'): [UID, Pre_cell, Post_cell], } ###


connection_types = {((None,'sim1'), (None,'sim1')): ['1490642434.41', 3, 5]}

pg.dbg()
app = pg.mkQApp()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
sweep_color = (0, 0, 0, 30)
avg_color = {'color': (255, 0, 255), 'width': 2}
holding_i = [-53, -60]
holding_e = [-65, -72]
holding = [-55, -72]
sign = '+'
scale_offset = (-10, -10)
scale_anchor = (0.45, 1)
sweep_threshold = 3

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
expt_cache = 'C:/Users/Stephanies/aisynphys/tools/expts_cache.pkl'
all_expts = cached_experiments()

test = PlotGrid()
test.set_shape(len(connection_types.keys()), 1)
test.show()
grid = PlotGrid()
if plot_trains is True:
    grid.set_shape(len(connection_types.keys()), 3)
    grid[0, 1].setTitle(title='50 Hz Train')
    grid[0, 2].setTitle(title='Exponential Deconvolution')
    tau = 15e-3
    lp = 1000
else:
    grid.set_shape(len(connection_types.keys()), 2)
grid.show()
row = 0
grid[0, 0].setTitle(title='First Pulse')

maxYpulse = []
maxYtrain = []
maxYdec = []
for row in range(len(connection_types)):
    connection_type = connection_types.keys()[row]
    if type(connection_type) is not tuple:
        connection_type = tuple(connection_type.split('-'))
    expt_id, pre_cell, post_cell = connection_types[connection_type] #all_connections[connection_type]
    expt = all_expts[expt_id]
    if expt.cells[pre_cell].cre_type in EXCITATORY_CRE_TYPES:
        holding = holding_e
        sign = '+'
    elif expt.cells[pre_cell].cre_type in INHIBITORY_CRE_TYPES:
        holding = holding_i
        sign = '-'
    pulse_response, artifact = get_response(expt, pre_cell, post_cell, analysis_type='pulse')
    sweep_list = response_filter(pulse_response, freq_range=[0, 50], holding_range=holding, pulse=True)
    n_sweeps = len(sweep_list[0])
    if n_sweeps > sweep_threshold:
        qc_list = pulse_qc(sweep_list, baseline=2, pulse=None, plot=pg.plot())
        qc_sweeps = len(qc_list)
        if qc_sweeps > sweep_threshold:
            avg_first_pulse = trace_avg(qc_list)
            avg_first_pulse.t0 = 0
            if plot_sweeps is True:
                for current_sweep in qc_list:
                    current_sweep.t0 = 0
                    bsub_trace = bsub(current_sweep)
                    trace_plot(bsub_trace, sweep_color, plot=grid[row, 0], x_range=[-2e-3, 27e-3])
            trace_plot(avg_first_pulse, avg_color, plot=grid[row, 0], x_range=[-2e-3, 27e-3])
            label = pg.LabelItem('%s, n = %d' % (connection_type, qc_sweeps))
            label.setParentItem(grid[row, 0].vb)
            label.setPos(50, 0)
            maxYpulse.append((row, grid[row,0].getAxis('left').range[1]))
            grid[row, 0].hideAxis('bottom')

        _, _, _, peak_t = get_amplitude(qc_list)
        all_amps = np.asarray(fail_rate(qc_list, sign=sign, peak_t=peak_t))
        # y = pg.pseudoScatter(all_amps, spacing=0.15)
        # test[row, 0].plot(all_amps, y, pen=None, symbol='o', symbolSize=8, symbolPen=(255, 255, 255, 200), symbolBrush=(0, 0, 255, 150))
        y,x = np.histogram(all_amps, bins=np.linspace(0, 2e-3, 40))
        test[row, 0].plot(x, y, stepMode=True, fillLevel=0, brush='k')
        test[row, 0].setLabels(bottom=('Vm', 'V'))
        test[row, 0].setXRange(0, 2e-3)
    else:
        print ("%s -> %s not enough sweeps for first pulse" % (connection_type[0], connection_type[1]))
    if row == len(connection_types) - 1:
        x_scale = pg.ScaleBar(size=10e-3, suffix='s')
        x_scale.setParentItem(grid[row, 0].vb)
        x_scale.anchor(scale_anchor, scale_anchor, offset=scale_offset)
    if plot_trains is True:
        train_responses, _ = get_response(expt, pre_cell, post_cell, analysis_type='train')
        train_sweep_list = response_filter(train_responses['responses'], freq_range=[50, 50], holding_range=holding, train=0)
        n_train_sweeps = len(train_sweep_list)
        if n_train_sweeps > sweep_threshold:
            dec_sweep_list = []
            for sweep in range(n_train_sweeps):
                train_sweep = train_sweep_list[sweep]
                train_base = bsub(train_sweep)
                dec_sweep = bessel_filter(exp_deconvolve(train_sweep, tau), lp)
                dec_base = bsub(dec_sweep)
                dec_sweep_list.append(dec_base)
                if plot_sweeps is True:
                    trace_plot(train_base, sweep_color, plot=grid[row, 1])
                    trace_plot(dec_base, sweep_color, plot=grid[row, 2])
            ind_avg = trace_avg(train_sweep_list)
            ind_avg.t0 = 0
            ind_dec = trace_avg(dec_sweep_list)
            ind_dec.t0 = 0
            trace_plot(ind_avg, avg_color, plot=grid[row, 1])
            trace_plot(ind_dec, avg_color, plot=grid[row, 2])
            label = pg.LabelItem('n = %d' % n_train_sweeps)
            label.setParentItem(grid[row, 1].vb)
            label.setPos(50, 0)
            grid[row, 1].label = label
            grid[row, 2].hideAxis('left')
            maxYtrain.append((row, grid[row, 1].getAxis('left').range[1]))
            maxYdec.append((row, grid[row, 2].getAxis('left').range[1]))
        else:
            print ("%s -> %s not enough sweeps for first pulse" % (connection_type[0], connection_type[1]))
            grid[row, 0].setYLink(grid[2, 1])
            grid[row, 1].hideAxis('bottom')
            grid[row, 2].hideAxis('bottom')
            grid[row, 1].hideAxis('left')
            grid[row, 2].hideAxis('left')
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

if sys.flags.interactive == 0:
    app.exec_()