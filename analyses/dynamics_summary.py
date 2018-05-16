from synapse_comparison import load_cache, summary_plot_pulse
from rep_connections import ee_connections, human_connections
from manuscript_figures import colors_mouse, feature_kw, colors_human
from scipy import stats
import numpy as np
import pyqtgraph as pg
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.synaptic_release import ReleaseModel

app = pg.mkQApp()
pg.dbg()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

add_model = False
grand_avg = False

cache_file = 'train_amps_human.pkl'
data = load_cache(cache_file)
expt_ids = load_cache('expt_ids_human.pkl')
model_amps = load_cache('model_amps.pkl')
feature_plt_ind = None
feature_plt_rec = None
symbols = ['d', 's', 'o', '+', 't']
color = colors_human

# model = ReleaseModel()
# model.Dynamics = {'Dep':1, 'Fac':0, 'UR':0, 'SMR':0, 'DSR':0}
# params = {(('2/3', 'unknown'), ('2/3', 'unknown')): [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#           ((None,'rorb'), (None,'rorb')): [0, 506.7, 0, 0, 0, 0, 0.22, 0, 0, 0, 0, 0],
#           ((None,'sim1'), (None,'sim1')): [0, 1213.8, 0, 0, 0, 0, 0.17, 0, 0, 0, 0, 0],
#           ((None,'tlx3'), (None,'tlx3')): [0, 319.4, 0, 0, 0, 0, 0.16, 0, 0, 0, 0, 0]}

ind = data[0]
rec = data[1]
freq = 10
delta =250
order = human_connections.keys()
delay_order = [250, 500, 1000, 2000, 4000]

ind_plt = PlotGrid()
ind_plt.set_shape(4, 1)
ind_plt.show()
rec_plt = PlotGrid()
rec_plt.set_shape(1, 2)
rec_plt.show()
ind_plt_scatter = pg.plot()
ind_plt_all = PlotGrid()
ind_plt_all.set_shape(1, 2)
ind_plt_all.show()
ind_50 = {}
ninth_pulse_250 = {}
gain_plot = pg.plot()


for t, type in enumerate(order):
    if type == (('2', 'unknown'), ('2', 'unknown')):
        continue
    if type == (('4', 'unknown'), ('4', 'unknown')):
        continue
    pulse_ratio = {}
    median_ratio = []
    sd_ratio = []
    gain = []
    for f, freqs in enumerate(ind[type]):
        # spike_times = (np.arange(12)*1/float(freqs[0])) * 1e3
        # spike_times[8:] += 250-(1/float(freqs[0]))*1e3
        # model_eval = model.eval(list(spike_times), params[type])
        if grand_avg is True:
            avg_pulses = np.mean(freqs[1], 0)
            avg_ratio = avg_pulses/avg_pulses[0]
            pulse_ratio[freqs[0]] = avg_ratio[7]
            gain.append((avg_ratio[8] - avg_ratio[7])/avg_ratio[8] *100)
        else:
            gain.append([(freqs[1][n, 8] - freqs[1][n, 7])/freqs[1][n, 8] * 100 for n in range(freqs[1].shape[0])])
            pulse_ratio[freqs[0]] = np.asarray([freqs[1][n, :]/freqs[1][n, 0] for n in range(freqs[1].shape[0])])
            avg_ratio = np.mean(pulse_ratio[freqs[0]], 0)
            sd_ratio = np.std(pulse_ratio[freqs[0]], 0)
            color2 = (color[t][0], color[t][1], color[t][2], 150)
            vals = np.hstack(pulse_ratio[freqs[0]][0])
            x = pg.pseudoScatter(vals, spacing=0.15)
            ind_plt_all[0, 1].plot(x, vals, pen=None, symbol=symbols[f], symbolSize=8, symbolPen=color[t],
                                   symbolBrush=color2)
            ind_plt_all[0, 1].setLabels(left=['8:1 Ratio', ''])
        ind_plt[t, 0].addLegend()
        ind_plt[t, 0].plot(avg_ratio, pen=color[t], symbol=symbols[f], symbolSize=10, symbolPen='k',
                        symbolBrush=color[t], name=('  %d Hz' % freqs[0]))
        if add_model is True:
            if type != (('2/3', 'unknown'), ('2/3', 'unknown')):
                model = model_amps[type][0][f]
                ind_plt[t, 0].plot(model, pen=color[t])
        ind_plt[t, 0].setXRange(0, 11)
        ind_plt[t, 0].setYRange(0, 1.5)
        ind_plt_all[0, 0].plot([f], [avg_ratio[7]], pen=None, symbol='o', symbolSize=15, symbolPen='k', symbolBrush=color[t])
        # err = pg.ErrorBarItem(x=np.asarray([f]), y=np.asarray([avg_ratio[7]]), height=np.asarray([sd_ratio[7]]), beam=0.1)
        # ind_plt_all[0, 0].addItem(err)

    ind_plt_all[0, 0].getAxis('bottom').setTicks([[(0, '10'), (1, '20'), (2, '50'), (3, '100')]])
    ind_plt_all[0, 0].setLabels(left=['8:1 Ratio', ''], bottom=['Frequency', 'Hz'])

    if grand_avg is False:
        feature_list = (pulse_ratio[10][:, 7], pulse_ratio[20][:, 7], pulse_ratio[50][:, 7], pulse_ratio[100][:, 7])
        labels = [['8:1 pulse ratio', ''], ['8:1 pulse ratio', ''], ['8:1 pulse ratio', ''], ['8:1 pulse ratio', '']]
        titles = ['10Hz', '20Hz', '50Hz', '100Hz']
        feature_plt_ind = summary_plot_pulse(feature_list, labels, titles, t, plot=feature_plt_ind,
                                         color=color[t], name=type)
        ind_50[type] = {}
        ind_50[type][50] = pulse_ratio[50][:, 7]
    gain_plot.plot([np.mean(n) for n in gain], pen=color[t], symbol='o', symbolSize=8, symbolPen='k',
                   symbolBrush=color[t])
    gain_plot.getAxis('bottom').setTicks([[(0, '10/4'), (1, '20/4'), (2, '50/4'), (3, '100/4')]])
    gain_plot.setLabels(left=['% change amplitude', ''], bottom=['Freq Change', ''])



    ninth_pulse = {}
    ind_pulses = []
    rec_pulse_ratio = {}
    rec_avg_ratio = []
    rec_perc = {}

    for d, delay in enumerate(rec[type]):
        recovery = []
        if grand_avg is True:
            avg_pulses = np.mean(delay[1], 0)
            rec_avg_ratio.append(avg_pulses/avg_pulses[0])
            nine = avg_pulses[8]/avg_pulses[0]
        else:
            rec_pulse_ratio[delay[0]] = np.asarray([delay[1][n, :]/delay[1][n, 0] for n in range(delay[1].shape[0])])
            rec_avg_ratio.append(np.mean(rec_pulse_ratio[delay[0]], 0))
            # rec_sd_ratio = np.std(rec_pulse_ratio[delay[0]], 0)
            nine = rec_pulse_ratio[delay[0]][:, 8]

        ninth_pulse[delay[0]] = nine

        if add_model is True:
            if type != (('2/3', 'unknown'), ('2/3', 'unknown')):
                model_rec = model_amps[type][1][d, 8:]
                rec_plt[0, 1].plot([d, d+0.2, d+0.4, d+0.6], model_rec, pen={'color': color[t], 'width': 2})

    grand_rec_ratio = np.mean(np.asarray(rec_avg_ratio), 0)
    rec_plt[0, 0].plot(grand_rec_ratio[:8]/grand_rec_ratio[0], pen=color[t], symbol='o',
                       symbolSize=10, symbolPen='k',symbolBrush=color[t])
    if add_model is True:
        if type != (('2/3', 'unknown'), ('2/3', 'unknown')):
            model_rec_ind = np.mean(model_amps[type][1], 0)[:8]
            rec_plt[0, 0].plot(model_rec_ind, pen={'color': color[t], 'width': 2})
    rec_plt[0, 0].getAxis('bottom').setTicks([[(0, '0'), (1, '20'), (2, '40'), (3, '60'), (4, '80'), (5, '100'), (6, '120'),
                                               (7, '140')]])
    rec_plt[0, 0].setYRange(0, 1.5)
    ninth_pulse_avg = [np.mean(ninth_pulse[delays]) for delays in delay_order]
    rec_plt[0, 1].plot(ninth_pulse_avg, pen=color[t], symbol='o', symbolSize=10, symbolPen='k',
                       symbolBrush=color[t])
    rec_plt[0, 1].setYRange(0, 1.5)
    rec_plt[0, 1].getAxis('bottom').setTicks([[(0, '250'), (1, '500'), (2, '1000'), (3, '2000'), (4, '4000')]])

    #
    if grand_avg is False:
        labels = labels = [['9:1 pulse ratio', ''], ['9:1 pulse ratio', ''], ['9:1 pulse ratio', ''],
                           ['9:1 pulse ratio', ''], ['9:1 pulse ratio', '']]
        titles = ['250ms', '500ms', '1000ms', '2000ms', '4000ms']
        feature_list = (ninth_pulse[250], ninth_pulse[500], ninth_pulse[1000], ninth_pulse[2000], ninth_pulse[4000])
        feature_plt_rec = summary_plot_pulse(feature_list, labels, titles, t, plot=feature_plt_rec,
                                         color=color[t], name=type)

    ninth_pulse_250[type] = {}
    ninth_pulse_250[type][250] = ninth_pulse[250]

if grand_avg is False:
    feature_kw(50, ind_50)
    feature_kw(250, ninth_pulse_250)