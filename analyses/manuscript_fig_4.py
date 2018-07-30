# *-* coding: utf-8 *-*
"""
2018 E-E manuscript fig 4:
Connectivity and its relationship to intersomatic distance
"""
from __future__ import print_function, division

from collections import OrderedDict

import numpy as np
import pyqtgraph as pg

from neuroanalysis.data import Trace
from multipatch_analysis.experiment_list import ExperimentList, cache_file


def write_csv(fh, data, description, units='connection probability %'):
    """Used to generate csv file accompanying figure.
    """
    if isinstance(data, Trace):
        write_csv(fh, data, description + "distance(um)")
        write_csv(fh, data, description + " %s" % units)
    else:
        cols = ['"' + description + '"'] + list(data)
        line = ','.join(map(str, cols))
        fh.write(line)
        fh.write('\n')


app = pg.mkQApp()

pg.dbg()
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

win = pg.GraphicsLayoutWidget()
win.show()
win.resize(1200, 600)



# set up distance plots
mouse_dist_plots = []
mouse_hist_plots = []
human_dist_plots = []
human_hist_plots = []
for row, plots in enumerate([(mouse_hist_plots, mouse_dist_plots), (human_hist_plots, human_dist_plots)]):
    hist_plots, dist_plots = plots
    xlabel = pg.LabelItem(u'distance (µm)')
    win.addItem(xlabel, row=row*3+2, col=0, colspan=5)
    for i in range(5):
        hist_plot = win.addPlot(row*3, i)
        dist_plot = win.addPlot(row*3+1, i)
        hist_plots.append(hist_plot)
        dist_plots.append(dist_plot)

        dist_plot.setXRange(20e-6, 180e-6)
        hist_plot.setXRange(20e-6, 180e-6)

        hist_plot.setMaximumHeight(40)
        dist_plot.setXLink(hist_plot)
        hist_plot.getAxis('bottom').hide()
        dist_plot.getAxis('bottom').setScale(1e6)
        dist_plot.getAxis('left').setScale(100)
        hist_plot.setXLink(mouse_hist_plots[0])

        if i == 0:
            dist_plot.setLabels(left='connection probability %')
            dist_plot.getAxis('left').setWidth(50)
            hist_plot.getAxis('left').setWidth(50)
        else:
            dist_plot.getAxis('left').setWidth(30)
            hist_plot.getAxis('left').setWidth(30)
            dist_plot.setYLink(dist_plots[0])
            hist_plot.setYLink(hist_plots[0])


app.processEvents()

all_expts = ExperimentList(cache=cache_file)


mouse_expts = all_expts.select(calcium='high', age='40-1000', organism='mouse')
human_expts = all_expts.select(organism='human')

mouse_ee_types = OrderedDict([
    (('2/3', 'unknown'), 'L23pyr'),
    ((None, 'rorb'), 'rorb'),
    ((None, 'tlx3'), 'tlx3'),
    ((None, 'sim1'), 'sim1'),
    ((None, 'ntsr1'), 'ntsr1'),
])

human_types = OrderedDict([
    (('1', 'unknown'), 'L1'),
    (('2', 'unknown'), 'L2'),
    (('3', 'unknown'), 'L3'), 
    (('4', 'unknown'), 'L4'), 
    (('5', 'unknown'), 'L5'), 
    (('6', 'unknown'), 'L6'),
])
human_types = OrderedDict([(typ, "L%s %s" % typ) for typ in human_types])

# mm = expts.matrix(mouse_ee_types, mouse_ee_types)
# hm = expts.matrix(human_types, human_types)


#populate mouse d/hist plots
#L2/3
plot, xvals, prop, upper, lower = mouse_expts.distance_plot([('2/3', None)], [('2/3', None)], plots=(mouse_dist_plots[0], None), color=(249, 144, 92))
connected, probed = mouse_expts.count_connections([('2/3', None)], [('2/3', None)])
bins = np.arange(0, 180e-6, 20e-6)
hist = np.histogram(probed, bins=bins)
mouse_hist_plots[0].plot(hist[1], hist[0], stepMode=True, fillLevel=0, brush=(249, 144, 92, 80))

fig4B = open('fig4B.csv', 'wb')

write_csv(fig4B, hist[1], "Figure4B, L2/3 histogram values")
write_csv(fig4B, hist[0], "Figure4B, L2/3 histogram bin edges")
write_csv(fig4B, xvals, "Figure 4B, L2/3 distance plot x vals")
write_csv(fig4B, prop, "Figure 4B, L2/3 distance plot trace")
write_csv(fig4B, upper, "Figure 4B, L2/3 distance plot upper CI")
write_csv(fig4B, lower, "Figure 4B, L2/3 distance plot x vals")

#rorb
mouse_expts.distance_plot('rorb', 'rorb', plots=(mouse_dist_plots[1], None), color=(100, 202, 103))
connected, probed = mouse_expts.count_connections('rorb', 'rorb')
bins = np.arange(0, 180e-6, 20e-6)
hist = np.histogram(probed, bins=bins)
mouse_hist_plots[1].plot(hist[1], hist[0], stepMode=True, fillLevel=0, brush=(100, 202, 103, 25))

write_csv(fig4B, hist[1], "Figure4B, rorb histogram values")
write_csv(fig4B, hist[0], "Figure4B, rorb histogram bin edges")
write_csv(fig4B, xvals, "Figure 4B, rorb distance plot x vals")
write_csv(fig4B, prop, "Figure 4B, rorb distance plot trace")
write_csv(fig4B, upper, "Figure 4B, rorb distance plot upper CI")
write_csv(fig4B, lower, "Figure 4B, rorb distance plot x vals")

#sim1
mouse_expts.distance_plot('sim1', 'sim1', plots=(mouse_dist_plots[2], None), color=(45, 77, 247))
connected, probed = mouse_expts.count_connections('sim1', 'sim1')
hist = np.histogram(probed, bins=bins)
mouse_hist_plots[2].plot(hist[1], hist[0], stepMode=True, fillLevel=0, brush=(45, 77, 247, 80))

write_csv(fig4B, hist[1], "Figure4B, sim1 histogram values")
write_csv(fig4B, hist[0], "Figure4B, sim1 histogram bin edges")
write_csv(fig4B, xvals, "Figure 4B, sim1 distance plot x vals")
write_csv(fig4B, prop, "Figure 4B, sim1 distance plot trace")
write_csv(fig4B, upper, "Figure 4B, sim1 distance plot upper CI")
write_csv(fig4B, lower, "Figure 4B, sim1 distance plot x vals")

#tlx3
mouse_expts.distance_plot('tlx3', 'tlx3', plots=(mouse_dist_plots[3], None), color=(81, 221, 209))
connected, probed = mouse_expts.count_connections('tlx3', 'tlx3')
hist = np.histogram(probed, bins=bins)
mouse_hist_plots[3].plot(hist[1], hist[0], stepMode=True, fillLevel=0, brush=(81, 221, 209, 80))

write_csv(fig4B, hist[1], "Figure4B, tlx3 histogram values")
write_csv(fig4B, hist[0], "Figure4B, tlx3 histogram bin edges")
write_csv(fig4B, xvals, "Figure 4B, tlx3 distance plot x vals")
write_csv(fig4B, prop, "Figure 4B, tlx3 distance plot trace")
write_csv(fig4B, upper, "Figure 4B, tlx3 distance plot upper CI")
write_csv(fig4B, lower, "Figure 4B, tlx3 distance plot x vals")

#ntsr1
mouse_expts.distance_plot('ntsr1', 'ntsr1', plots=(mouse_dist_plots[4], None), color=(153, 51, 255))
connected, probed = mouse_expts.count_connections('ntsr1', 'ntsr1')
hist = np.histogram(probed, bins=bins)
mouse_hist_plots[4].plot(hist[1], hist[0], stepMode=True, fillLevel=0, brush=(153, 51, 255, 80))

write_csv(fig4B, hist[1], "Figure4B, ntsr1 histogram values")
write_csv(fig4B, hist[0], "Figure4B, ntsr1 histogram bin edges")
write_csv(fig4B, xvals, "Figure 4B, ntsr1 distance plot x vals")
write_csv(fig4B, prop, "Figure 4B, ntsr1 distance plot trace")
write_csv(fig4B, upper, "Figure 4B, ntsr1 distance plot upper CI")
write_csv(fig4B, lower, "Figure 4B, ntsr1 distance plot x vals")

fig4B.close()

#populate human d/hist plots
#L2

fig4D = open('fig4D.csv', 'wb')

human_expts.distance_plot([('2', None)], [('2', None)], plots=(human_dist_plots[0], None), color=(247, 118, 118))
connected, probed = human_expts.count_connections([('2', None)], [('2', None)])
bins = np.arange(0, 180e-6, 20e-6)
hist = np.histogram(probed, bins=bins)
human_hist_plots[0].plot(hist[1], hist[0], stepMode=True, fillLevel=0, brush=(247, 118, 118, 80))

write_csv(fig4D, hist[1], "Figure4D, L2 histogram values")
write_csv(fig4D, hist[0], "Figure4D, L2 histogram bin edges")
write_csv(fig4D, xvals, "Figure 4D, L2 distance plot x vals")
write_csv(fig4D, prop, "Figure 4D, L2 distance plot trace")
write_csv(fig4D, upper, "Figure 4D, L2 distance plot upper CI")
write_csv(fig4D, lower, "Figure 4D, L2 distance plot x vals")

#L3
human_expts.distance_plot([('3', None)], [('3', None)], plots=(human_dist_plots[1], None), color=(246, 197, 97))
connected, probed = human_expts.count_connections([('3', None)], [('3', None)])
bins = np.arange(0, 180e-6, 20e-6)
hist = np.histogram(probed, bins=bins)
human_hist_plots[1].plot(hist[1], hist[0], stepMode=True, fillLevel=0, brush=(246, 197, 97, 80))

write_csv(fig4D, hist[1], "Figure4D, L3 histogram values")
write_csv(fig4D, hist[0], "Figure4D, L3 histogram bin edges")
write_csv(fig4D, xvals, "Figure 4D, L3 distance plot x vals")
write_csv(fig4D, prop, "Figure 4D, L3 distance plot trace")
write_csv(fig4D, upper, "Figure 4D, L3 distance plot upper CI")
write_csv(fig4D, lower, "Figure 4D, L3 distance plot x vals")

#L4
human_expts.distance_plot([('4', None)], [('4', None)], plots=(human_dist_plots[2], None), color=(100, 202, 103))
connected, probed = human_expts.count_connections([('4', None)], [('4', None)])
bins = np.arange(0, 180e-6, 20e-6)
hist = np.histogram(probed, bins=bins)
human_hist_plots[2].plot(hist[1], hist[0], stepMode=True, fillLevel=0, brush=(100, 202, 103, 80))

write_csv(fig4D, hist[1], "Figure4D, L4 histogram values")
write_csv(fig4D, hist[0], "Figure4D, L4 histogram bin edges")
write_csv(fig4D, xvals, "Figure 4D, L4 distance plot x vals")
write_csv(fig4D, prop, "Figure 4D, L4 distance plot trace")
write_csv(fig4D, upper, "Figure 4D, L4 distance plot upper CI")
write_csv(fig4D, lower, "Figure 4D, L4 distance plot x vals")

#L5
human_expts.distance_plot([('5', None)], [('5', None)], plots=(human_dist_plots[3], None), color=(107, 155, 250))
connected, probed = human_expts.count_connections([('5', None)], [('5', None)])
bins = np.arange(0, 180e-6, 20e-6)
hist = np.histogram(probed, bins=bins)
human_hist_plots[3].plot(hist[1], hist[0], stepMode=True, fillLevel=0, brush=(107, 155, 250, 80))

write_csv(fig4D, hist[1], "Figure4D, L5 histogram values")
write_csv(fig4D, hist[0], "Figure4D, L5 histogram bin edges")
write_csv(fig4D, xvals, "Figure 4D, L5 distance plot x vals")
write_csv(fig4D, prop, "Figure 4D, L5 distance plot trace")
write_csv(fig4D, upper, "Figure 4D, L5 distance plot upper CI")
write_csv(fig4D, lower, "Figure 4D, L5 distance plot x vals")

fig4D.close()

#L6
human_expts.distance_plot([('6', None)], [('6', None)], plots=(human_dist_plots[4], None), color=(153, 51, 255))
connected, probed = human_expts.count_connections([('6', None)], [('6', None)])
bins = np.arange(0, 180e-6, 20e-6)
hist = np.histogram(probed, bins=bins)
human_hist_plots[4].plot(hist[1], hist[0], stepMode=True, fillLevel=0, brush=(153, 51, 255, 80))

#fix plot sizing
mouse_dist_plots[0].setXRange(0,180e-6)
mouse_dist_plots[0].setYRange(0,0.3)
mouse_hist_plots[0].setYRange(0,215)

human_dist_plots[0].setXRange(0,180e-6)
human_dist_plots[0].setYRange(0,0.6)
human_hist_plots[0].setYRange(0,35)

