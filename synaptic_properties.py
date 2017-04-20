import pyqtgraph as pg
from experiment_list import ExperimentList
from connection_detection import MultiPatchSyncRecAnalyzer, EvokedResponseGroup, trace_mean
from neuroanalysis.ui.plot_grid import PlotGrid

app = pg.mkQApp()
pg.dbg()

all_expts = ExperimentList(cache='expts_cache.pkl')

def plot_trace_average(expts, pre_type, post_type, avg_plot, ind_plot=None, clamp_mode='ic', stim_filter='50Hz'):
    all_responses = []
    n_traces = 0
    
    for conn in all_expts.connection_summary()[::-1]:
        pre_cell, post_cell = conn['cells']
        if pre_cell.cre_type != pre_type or post_cell.cre_type != post_type:
            continue
        expt = conn['expt']
        pre_id, post_id = pre_cell.cell_id-1, post_cell.cell_id-1
        print expt, pre_id, post_id
        #with expt.data as data:
        data = expt.data
    
        responses = EvokedResponseGroup(pre_id, post_id)
        for sweep in data.contents:
            devs = sweep.devices
            if pre_id not in devs or post_id not in devs:
                continue
            pre_rec = sweep[pre_id]
            if stim_filter not in pre_rec.meta['stim_name']:
                continue
            post_rec = sweep[post_id]
            if post_rec.clamp_mode != clamp_mode:
                continue
            
            sa = MultiPatchSyncRecAnalyzer.get(sweep)
            resp = sa.get_pulse_response(pre_rec, post_rec, 0, 7)
            base = None
            if resp is None or resp.duration < 200e-3:
                continue
            responses.add(resp, base)
        
        print "   got %d sweeps" % len(responses)
        
        if len(responses) < 25:
            continue
        
        avg = responses.mean()
        n_traces += len(responses)
        all_responses.append(avg)
        if ind_plot is not None:
            ind_plot.plot(avg.time_values, avg.data)
            ind_plot.setTitle('%s - %s  (%d synapses, %d traces)' % (pre_type, post_type, len(all_responses), n_traces))

        app.processEvents()

    avg_plot.setTitle('%s - %s  (%d synapses, %d traces)' % (pre_type, post_type, len(all_responses), n_traces))
    if len(all_responses) > 0:
        avg = trace_mean(all_responses)
        avg_plot.plot(avg.time_values, avg.data, antialias=True)
        avg_plot.setLabels(left=('Average response', 'V'), bottom=('Time', 's'))


if __name__ == '__main__':
    #avg_plot = pg.plot()
    #ind_plot = pg.plot()
    #plot_trace_average(all_expts, 'tlx3', 'tlx3', avg_plot, ind_plot, clamp_mode='ic', stim_filter='50Hz')
    

    types = ['tlx3', 'sim1', 'pvalb', 'sst', 'vip']
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
            plot_trace_average(all_expts, pre_type, post_type, avg_plot, ind_plot, clamp_mode='ic', stim_filter='50Hz')
