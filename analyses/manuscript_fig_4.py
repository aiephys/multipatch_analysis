# *-* coding: utf-8 *-*
"""
2018 E-E manuscript fig 4:
Connectivity and its relationship to intersomatic distance
"""
from __future__ import print_function, division

from collections import OrderedDict

import numpy as np
import pyqtgraph as pg
import scipy.stats

from neuroanalysis.data import TSeries
from aisynphys.experiment_list import ExperimentList, cache_file
from aisynphys.cell_class import CellClass, classify_cells
from aisynphys.connectivity import query_pairs, measure_connectivity
from aisynphys.database import database as db
from aisynphys.ui.graphics import distance_plot


def write_csv(fh, data, description, units='connection probability %'):
    """Used to generate csv file accompanying figure.
    """
    if isinstance(data, TSeries):
        write_csv(fh, data, description + "distance(um)")
        write_csv(fh, data, description + " %s" % units)
    else:
        cols = ['"' + description + '"'] + list(data)
        line = ','.join(map(str, cols))
        fh.write(line)
        fh.write('\n')


# Set up UI
app = pg.mkQApp()
# pg.dbg()

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

win = pg.GraphicsLayoutWidget()
win.show()
win.resize(800, 480)

# set up connectivity plots
mouse_conn_plot = win.addPlot(0, 0, rowspan=3, labels={'left': 'connection probability %', 'bottom': ''})
human_conn_plot = win.addPlot(3, 0, rowspan=3, labels={'left': 'connection probability %', 'bottom': ''})
human_conn_plot.setXLink(mouse_conn_plot)
mouse_conn_plot.setXRange(-0.5, 4.5, padding=0)
for plt in mouse_conn_plot, human_conn_plot:
    plt.setFixedWidth(220)
    plt.setYRange(0, 0.3)
    plt.getAxis('left').setStyle(tickLength=0)
    plt.getAxis('bottom').setStyle(tickLength=0)
    plt.getAxis('left').setScale(100)


# set up distance plots
mouse_dist_plots = []
mouse_hist_plots = []
human_dist_plots = []
human_hist_plots = []
for row, plots in enumerate([(mouse_hist_plots, mouse_dist_plots), (human_hist_plots, human_dist_plots)]):
    hist_plots, dist_plots = plots
    xlabel = pg.LabelItem(u'distance (µm)')
    xlabel.setFixedHeight(20)
    win.addItem(xlabel, row=row*3+2, col=1, colspan=5)
    for i in range(5):
        hist_plot = win.addPlot(row*3, i+1)
        dist_plot = win.addPlot(row*3+1, i+1)
        hist_plots.append(hist_plot)
        dist_plots.append(dist_plot)

        hist_plot.setMaximumHeight(40)
        dist_plot.setXLink(hist_plot)
        hist_plot.getAxis('bottom').hide()
        dist_plot.getAxis('bottom').setScale(1e6)
        dist_plot.getAxis('left').setScale(100)
        hist_plot.setXLink(mouse_hist_plots[0])

        hist_plot.hideAxis('bottom')

        if i == 0:
            dist_plot.setLabels(left='connection probability %')
            dist_plot.getAxis('left').setWidth(50)
            hist_plot.getAxis('left').setWidth(50)
        else:
            dist_plot.getAxis('left').setWidth(0)
            hist_plot.getAxis('left').setWidth(0)
            dist_plot.setYLink(dist_plots[0])
            hist_plot.setYLink(hist_plots[0])

            dist_plot.getAxis('left').setStyle(showValues=False, tickTextWidth=0, tickTextHeight=0, tickLength=0)
            hist_plot.getAxis('left').setStyle(showValues=False, tickTextWidth=0, tickTextHeight=0, tickLength=0)
        dist_plot.getAxis('left').setStyle(tickLength=0)
        hist_plot.getAxis('left').setStyle(tickLength=0)
        dist_plot.getAxis('bottom').setStyle(tickLength=0)

        dist_plot.vb.setMinimumWidth(80)
        hist_plot.vb.setMinimumWidth(80)

mouse_hist_plots[0].getAxis('left').setTicks([[(0, '0'), (200, '200')]])
human_hist_plots[0].getAxis('left').setTicks([[(0, '0'), (50, '50')]])


app.processEvents()

mouse_dist_plots[0].setXRange(20e-6, 180e-6, padding=0)
mouse_dist_plots[0].setYRange(0, 0.4, padding=0)
human_dist_plots[0].setYRange(0, 0.4, padding=0)
mouse_hist_plots[0].setYRange(0, 200)
human_hist_plots[0].setYRange(0, 50)


app.processEvents()


# define CellClass, short name, and color for each class
mouse_classes = [
    (CellClass(target_layer='2/3', pyramidal=True), 'L2/3', (249, 144, 92)),
    (CellClass(cre_type='rorb'), 'Rorb', (100, 202, 103)),
    (CellClass(cre_type='tlx3'), 'Tlx3', (81, 221, 209)),
    (CellClass(cre_type='sim1'), 'Sim1', (45, 77, 247)),
    (CellClass(cre_type='ntsr1'), 'Ntsr1', (153, 51, 255)),
]

human_classes = [
    (CellClass(target_layer='2', pyramidal=True), 'L2', (247, 118, 118)),
    (CellClass(target_layer='3', pyramidal=True), 'L3', (246, 197, 97)),
    (CellClass(target_layer='4', pyramidal=True), 'L4', (100, 202, 103)),
    (CellClass(target_layer='5', pyramidal=True), 'L5', (107, 155, 250)),
    (CellClass(target_layer='6', pyramidal=True), 'L6', (153, 51, 255)),
]


session = db.session()


max_distance = None#100e-6

# analyze connectivity <100 um for mouse
mouse_pairs = query_pairs(acsf="2mM Ca & Mg", age=(40, None), species='mouse', distance=(None, max_distance), session=session).all()
mouse_groups = classify_cells([c[0] for c in mouse_classes], pairs=mouse_pairs)
mouse_results = measure_connectivity(mouse_pairs, mouse_groups)

# analyze connectivity < 100 um for human
human_pairs = query_pairs(species='human', distance=(None, max_distance), session=session).all()
human_groups = classify_cells([c[0] for c in human_classes], pairs=human_pairs)
human_results = measure_connectivity(human_pairs, human_groups)

csv_file = open("manuscript_fig_4.csv", 'wb')

# plot connectivity < 100 um
for conn_plot, results, cell_classes, fig_letter in [
        (mouse_conn_plot, mouse_results, mouse_classes, 'A'),
        (human_conn_plot, human_results, human_classes, 'C')]:

    # set axis labels
    class_names = [class_info[1] for class_info in cell_classes]
    conn_plot.getAxis('bottom').setTicks([list(enumerate(class_names))])

    # make arrays of results for plotting
    conn_data = np.empty((len(cell_classes), 3))
    n_probed = []
    n_connected = []
    brushes = []
    pens = []
    for i,class_info in enumerate(cell_classes):
        cell_class, _, color = class_info
        class_results = results[(cell_class, cell_class)]
        conn_data[i] = class_results['connection_probability']
        n_probed.append(class_results['n_probed'])
        n_connected.append(class_results['n_connected'])
        print("Cell class: %s  connected: %d  probed: %d  probability: %0.3f  min_ci: %0.3f  max_ci: %0.3f" % (
            cell_class.name, class_results['n_connected'], class_results['n_probed'], 
            class_results['connection_probability'][0], class_results['connection_probability'][1], class_results['connection_probability'][2],
        ))
        brushes.append(pg.mkBrush(color))

        # Add confidence interval line for this class
        errbar = pg.QtGui.QGraphicsLineItem(i, class_results['connection_probability'][1], i, class_results['connection_probability'][2])
        errbar.setPen(pg.mkPen(color, width=3))
        conn_plot.addItem(errbar)

    write_csv(csv_file, class_names, "Figure 4%s cell classes" % fig_letter)
    write_csv(csv_file, n_connected, "Figure 4%s n connected pairs < 100um" % fig_letter)
    write_csv(csv_file, n_probed, "Figure 4%s n probed pairs < 100um" % fig_letter)

    # plot connection probability points
    conn_plot.plot(conn_data[:,0], pen=None, symbol='o', symbolBrush=brushes, symbolPen='k')
    write_csv(csv_file, conn_data[:, 0], "Figure 4%s connection probability < 100um" % fig_letter)

    # write confidence intervals into csv
    write_csv(csv_file, conn_data[:, 1], "Figure 4%s connection probability < 100um lower 95%% confidence interval" % fig_letter)
    write_csv(csv_file, conn_data[:, 2], "Figure 4%s connection probability < 100um upper 95%% confidence interval" % fig_letter)


# statistical tests
# Is mouse L6 really different from other layers?
compare_classes = ['L2/3 pyr', 'rorb', 'tlx3', 'sim1']
print(u"========================\nMouse L6 connectivity < 100 µm vs %d classes (bonferroni corrected p=%0.3f:" % (len(compare_classes), 0.05 / len(compare_classes)))
a, b = mouse_results[('ntsr1', 'ntsr1')]['n_connected'], mouse_results[('ntsr1', 'ntsr1')]['n_probed']
for comp_class in compare_classes:
    c, d = mouse_results[(comp_class, comp_class)]['n_connected'], mouse_results[(comp_class, comp_class)]['n_probed']
    table = [[a, b], [c, d]]
    fe_r, fe_p = scipy.stats.fisher_exact(table)
    print("   {:10s} {:20s} p={:g}".format(comp_class, str(table), fe_p))
        
# Is human L4 really different from other layers?
compare_classes = ['L2 pyr', 'L3 pyr', 'L5 pyr']
print(u"========================\nHuman L4 connectivity < 100 µm vs %d classes (bonferroni corrected p=%0.3f:" % (len(compare_classes), 0.05 / len(compare_classes)))
a, b = human_results[('L4 pyr', 'L4 pyr')]['n_connected'], human_results[('L4 pyr', 'L4 pyr')]['n_probed']
for comp_class in compare_classes:
    c, d = human_results[(comp_class, comp_class)]['n_connected'], human_results[(comp_class, comp_class)]['n_probed']
    table = [[a, b], [c, d]]
    fe_r, fe_p = scipy.stats.fisher_exact(table)
    print("   {:10s} {:20s} p={:g}".format(comp_class, str(table), fe_p))
        

# does mouse connectivity differ from human?
# (this is a weird analysis given that we have sampled the classes unevenly)
table = [[0, 0], [0, 0]]
for i, species_opts in enumerate([(mouse_classes, mouse_results), (human_classes, human_results)]):
    cell_classes, results = species_opts
    for class_opts in cell_classes:
        cell_class = class_opts[0]
        a, b = results[(cell_class, cell_class)]['n_connected'], results[(cell_class, cell_class)]['n_probed']
        table[i][0] += a
        table[i][1] += b
print(u"========================\nHuman vs mouse overall connectivity < 100 µm:")
fe_r, fe_p = scipy.stats.fisher_exact(table)
print("    {:s} r={:0.2f} p={:g}".format(str(table), fe_r, fe_p))




# analyze connectivity (all distances) for mouse
#   note: the distance filter removes pairs that don't have a distance reported
mouse_pairs = query_pairs(acsf="2mM Ca & Mg", age=(40, None), distance=(None, float('inf')), species='mouse', session=session).all()
mouse_groups = classify_cells([c[0] for c in mouse_classes], pairs=mouse_pairs)
mouse_results = measure_connectivity(mouse_pairs, mouse_groups)

# analyze connectivity (all distances) for human
human_pairs = query_pairs(species='human', distance=(None, float('inf')), session=session).all()
human_groups = classify_cells([c[0] for c in human_classes], pairs=human_pairs)
human_results = measure_connectivity(human_pairs, human_groups)

# plot connectivity vs distance
for dist_plots, hist_plots, results, cell_classes, fig_letter in [
        (mouse_dist_plots, mouse_hist_plots, mouse_results, mouse_classes, 'B'),
        (human_dist_plots, human_hist_plots, human_results, human_classes, 'D')]:

    # iterate over cell classes
    for i, class_info in enumerate(cell_classes):
        cell_class, class_name, color = class_info

        # Add text label
        label = pg.LabelItem(class_name)
        label.setFixedWidth(80)
        label.setParentItem(dist_plots[i].vb)

        class_results = results[(cell_class, cell_class)]

        probed_pairs = class_results['probed_pairs']
        connected_pairs = class_results['connected_pairs']

        probed_distance = [p.distance for p in probed_pairs]
        connections = [(p in connected_pairs) for p in probed_pairs]

        plot, xvals, prop, upper, lower = distance_plot(connections, probed_distance, plots=(dist_plots[i], None), color=color, window=40e-6, spacing=40e-6)

        bins = np.arange(0, 180e-6, 20e-6)
        hist = np.histogram(probed_distance, bins=bins)
        hist_plots[i].plot(hist[1], hist[0], stepMode=True, fillLevel=0, brush=color + (80,))

        write_csv(csv_file, hist[1], "Figure 4%s, %s histogram values" % (fig_letter, class_name))
        write_csv(csv_file, hist[0], "Figure 4%s, %s histogram bin edges" % (fig_letter, class_name))
        write_csv(csv_file, xvals, "Figure 4%s, %s distance plot x vals" % (fig_letter, class_name))
        write_csv(csv_file, prop, "Figure 4%s, %s distance plot trace" % (fig_letter, class_name))
        write_csv(csv_file, upper, "Figure 4%s, %s distance plot upper CI" % (fig_letter, class_name))
        write_csv(csv_file, lower, "Figure 4%s, %s distance plot x vals" % (fig_letter, class_name))


