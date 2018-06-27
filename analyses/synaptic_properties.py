import pyqtgraph as pg
from multipatch_analysis.experiment_list import cached_experiments
from multipatch_analysis.connection_detection import MultiPatchExperimentAnalyzer, MultiPatchSyncRecAnalyzer, EvokedResponseGroup, trace_mean
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import Trace


def find_connections(expts, pre_type, post_type):
    """Iterate over all connections having a certain pre- and post-cre type.
    Yields, expt, pre_is, post_id for each connection.
    """
    for conn in all_expts.connection_summary()[::-1]:
        pre_cell, post_cell = conn['cells']
        if pre_cell.cre_type != pre_type or post_cell.cre_type != post_type:
            continue
        expt = conn['expt']
        pre_id, post_id = pre_cell.cell_id-1, post_cell.cell_id-1
        yield expt, pre_id, post_id


def find_sweeps(expt, pre_id, post_id, clamp_mode=None, stim_filter=None):
    """Iterate over all sweeps in an experiment having a particular pre- and post-
    device ID
    """
    for sweep in expt.contents:
        devs = sweep.devices
        if pre_id not in devs or post_id not in devs:
            continue
        pre_rec = sweep[pre_id]
        if stim_filter is not None and stim_filter not in pre_rec.stimulus.description:
            continue
        post_rec = sweep[post_id]
        if clamp_mode is not None and post_rec.clamp_mode != clamp_mode:
            continue
        yield sweep, pre_rec, post_rec



def plot_trace_average(expts, pre_type, post_type, avg_plot, ind_plot=None, clamp_mode='ic', stim_filter='50Hz', min_traces=25):
    all_responses = []
    n_traces = 0
    
    for expt, pre_id, post_id in find_connections(expts, pre_type, post_type):
        print expt, pre_id, post_id
        #with expt.data as data:
        data = expt.data
    
        responses = EvokedResponseGroup(pre_id, post_id)
        for sweep, pre_rec, post_rec in find_sweeps(data, pre_id, post_id, clamp_mode, stim_filter):
            sa = MultiPatchSyncRecAnalyzer.get(sweep)
            resp = sa.get_pulse_response(pre_rec, post_rec, 0, 7)
            base = None
            if resp is None or resp.duration < 200e-3:
                continue
            responses.add(resp, base)
        
        print "   got %d sweeps" % len(responses)
        
        if len(responses) < min_traces:
            continue
        
        avg = responses.mean()
        avg.t0 = 0
        n_traces += len(responses)
        all_responses.append(avg)
        if ind_plot is not None:
            ind_plot.plot(avg.time_values, avg.data)
            ind_plot.setTitle('%s - %s  (%d synapses, %d traces)' % (pre_type, post_type, len(all_responses), n_traces))

        pg.QtGui.QApplication.processEvents()

    avg_plot.setTitle('%s - %s  (%d synapses, %d traces)' % (pre_type, post_type, len(all_responses), n_traces))
    if len(all_responses) > 0:
        avg = trace_mean(all_responses)
        avg_plot.plot(avg.time_values, avg.data, antialias=True)
        avg_plot.setLabels(left=('Average response', 'V'), bottom=('Time', 's'))


def trace_average_matrix(expts, **kwds):
    types = ['tlx3', 'sim1', 'pvalb', 'sst', 'vip']
    results = {}
    plots = PlotGrid()
    plots.set_shape(len(types), len(types))
    indplots = PlotGrid()
    indplots.set_shape(len(types), len(types))
    plots.show()
    indplots.show()
    for i, pre_type in enumerate(types):
        for j, post_type in enumerate(types):
            avg_plot = plots[i, j]
            ind_plot = indplots[i, j]
            avg = plot_trace_average(all_expts, pre_type, post_type, avg_plot, ind_plot, **kwds)
            results[(pre_type, post_type)] = avg
    return results


if __name__ == '__main__':
    #avg_plot = pg.plot()
    #ind_plot = pg.plot()
    #plot_trace_average(all_expts, 'tlx3', 'tlx3', avg_plot, ind_plot, clamp_mode='ic', stim_filter='50Hz')
    app = pg.mkQApp()
    pg.dbg()

    all_expts = cached_experiments()
    
    results = trace_average_matrix(all_expts, clamp_mode='ic', stim_filter='50Hz', min_traces=25)
