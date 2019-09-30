from synapse_comparison import load_cache
import pyqtgraph as pg
from scipy import stats
from aisynphys.constants import EXCITATORY_CRE_TYPES, INHIBITORY_CRE_TYPES
from manuscript_figures import get_color
import numpy as np
import operator

pg.dbg()
app = pg.mkQApp()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

feature_cache = 'feature_cache_file.pkl'
features = load_cache(feature_cache)
freqs = [10, 20, 50, 100, 200]
rec_t = [250, 500, 1000, 2000, 4000]
symbols = ['o','s','t','d','+', 'o', 's','t']
connections = features['Amplitudes'].keys()

avg_amp = {}
abs_avg_amp = {}
for key, amp in features['Amplitudes'].items():
    abs_avg_amp[key] = abs(np.mean(amp))
    avg_amp[key] = np.nanmean(amp)

abs_sort_amp = sorted(abs_avg_amp.items(), key = operator.itemgetter(1), reverse=True)

i = 0
amp_plot = pg.plot()
amp_plot2 = pg.plot()
label = []
avg_amps = [a[1] for a in abs_sort_amp]
amps_sem = []
bar = pg.BarGraphItem(x=range(len(avg_amps)), height=avg_amps, width=0.7)
amp_plot.addItem(bar)
for key, avg_amp in abs_sort_amp:
    color= get_color (key[0], key[1])
    connection = ('%s->%s' % (key[0], key[1]))
    label.append((i, connection))
    amps = [abs(a) for a in features['Amplitudes'][key]]
    amps_sem.append(stats.sem(amps))
    dx = pg.pseudoScatter(np.array(amps).astype(float), 0.3, bidir=True)
    amp_plot.plot((0.3 * dx / dx.max()) + i, amps, pen=None, symbol='o', symbolSize=8, symbolBrush=color, symbolPen='w')
    amp_plot2.plot((0.3 * dx / dx.max()) + i, amps, pen=None, symbol='x', symbolSize=8, symbolBrush=color,
                  symbolPen=None)
    amp_plot2.plot([i], [avg_amp], pen=None, symbol='o', symbolBrush=color, symbolPen='w', symbolSize=15)
    i += 1

amp_plot.setLabels(left=('Abs Amplitude', 'V'))
amp_plot.getAxis('bottom').setTicks([label])
amp_plot.setYRange(0, 3e-3)
err = pg.ErrorBarItem(x=np.array(range(len(amps_sem))), y=np.array(avg_amps), height=np.array(amps_sem), beam=0.3)
amp_plot.addItem(err)

amp_plot2.setLabels(left=('Abs Amplitude', 'V'))
amp_plot2.getAxis('bottom').setTicks([label])
amp_plot2.setYRange(0, 3e-3)
err2 = pg.ErrorBarItem(x=np.array(range(len(amps_sem))), y=np.array(avg_amps), height=np.array(amps_sem), beam=0.3)
amp_plot2.addItem(err2)


ind_index_avg = np.zeros((20, 5))
ind_plot = pg.plot()
for f, freq in enumerate(features['Induction'].keys()):
    for c, type in enumerate(connections):
        if type in features['Induction'][freq].keys():
            ind_index_avg[c, f] = np.mean(features['Induction'][freq][type])

ind_sort = -np.sort(-ind_index_avg[:, 1])

rec_index_avg = np.zeros((20, 5))
for t, delta in enumerate(features['Recovery'].keys()):
    for c, type in enumerate(connections):
        if type in features['Recovery'][delta].keys():
            rec_index_avg[c, t] = np.mean(features['Recovery'][delta][type])



# ee_ind = pg.plot()
# ee_ind.addLegend()
# ee_ind.getAxis('bottom').setLogMode(False)
# ei_ind = pg.plot()
# ei_ind.addLegend()
# ei_ind.getAxis('bottom').setLogMode(False)
# ie_ind = pg.plot()
# ie_ind.addLegend()
# ie_ind.getAxis('bottom').setLogMode(False)
# ii_ind = pg.plot()
# ii_ind.addLegend()
# ii_ind.getAxis('bottom').setLogMode(False)
#
# ee_rec = pg.plot()
# ee_rec.addLegend()
# ee_rec.getAxis('bottom').setLogMode(False)
# ei_rec = pg.plot()
# ei_rec.addLegend()
# ei_rec.getAxis('bottom').setLogMode(False)
# ie_rec = pg.plot()
# ie_rec.addLegend()
# ie_rec.getAxis('bottom').setLogMode(False)
# ii_rec = pg.plot()
# ii_rec.addLegend()
# ii_rec.getAxis('bottom').setLogMode(False)
# for c, type in enumerate(connections):
#     pre_type = type[0]
#     post_type = type[1]
#     color = get_color(pre_type, post_type)
#     if pre_type in EXCITATORY_CRE_TYPES and post_type in EXCITATORY_CRE_TYPES:
#         ind_plt = ee_ind
#         rec_plt = ee_rec
#     elif pre_type in EXCITATORY_CRE_TYPES and post_type in INHIBITORY_CRE_TYPES:
#         ind_plt = ei_ind
#         rec_plt = ei_rec
#     elif pre_type in INHIBITORY_CRE_TYPES and post_type in EXCITATORY_CRE_TYPES:
#         ind_plt = ie_ind
#         rec_plt = ie_rec
#     elif pre_type in INHIBITORY_CRE_TYPES and post_type in INHIBITORY_CRE_TYPES:
#         ind_plt = ii_ind
#         rec_plt = ii_rec
#     ind_plt.plot(freqs, ind_index_avg[c,], pen=color, symbol='o', symbolBrush=color, symbolSize=8, symbolPen=None,
#                  name=('%s -> %s' % (type[0], type[1])))
#     ind_plt.setLabels(left=('Ratio 8th/1st', ''), bottom=('Induction Frequency', 'Hz'))
#     rec_plt.plot(rec_t, rec_index_avg[c,], pen=color, symbol='o', symbolBrush=color, symbolSize=8, symbolPen=None,
#                  name=('%s -> %s' % (type[0], type[1])))
#     rec_plt.setLabels(left=('Ratio 9th/8th', ''), bottom=('Recovery Delta T', 's'))