from collections import OrderedDict
import pyqtgraph as pg
import numpy as np
import os, sys, pickle, tempfile, resource
from neuroanalysis.ui.plot_grid import PlotGrid
from multipatch_analysis.experiment_list import cached_experiments()
from multipatch_analysis.connection_detection import MultiPatchExperimentAnalyzer, EvokedResponseGroup, fit_psp
from synaptic_properties import find_connections


def get_pulse_responses(expts, pre_type, post_type, **kwds):
    all_responses = EvokedResponseGroup()
    for expt, pre_id, post_id in find_connections(expts, pre_type, post_type):
        print expt, pre_id, post_id
        #with expt.data as data:
        data = expt.data
    
        xa = MultiPatchExperimentAnalyzer.get(data)
        responses = xa.get_evoked_responses(pre_id, post_id, **kwds)
        avg = responses.bsub_mean()
        if avg is not None:
            all_responses.add(avg, None)
        pg.QtGui.QApplication.processEvents()
    
    return all_responses
        
        
def plot_pulse_average(expts, pre_type, post_type, avg_plot, ind_plot, **kwds):
    all_avgs = EvokedResponseGroup()
    ind_plot.setLabels(left=('%s-%s'%(pre_type, post_type), 'V'), bottom=('Time', 's'))
    avg_plot.setLabels(left=('%s-%s'%(pre_type, post_type), 'V'), bottom=('Time', 's'))
    responses = get_pulse_responses(expts, pre_type, post_type, **kwds)
    for resp in responses.responses:
        ind_plot.plot(resp.time_values, resp.data)
        all_avgs.add(Trace(resp.data, dt=resp.dt), None)
        
    avg = all_avgs.mean()
    if avg is not None:
        avg_plot.plot(avg.time_values, avg.data)
    return avg


def pulse_average_matrix(expts, **kwds):
    results = {}
    types = ['sim1', 'tlx3', 'pvalb', 'sst', 'vip']
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
            avg = plot_pulse_average(all_expts, pre_type, post_type, avg_plot, ind_plot, **kwds)
            pg.QtGui.QApplication.processEvents()
            results[(pre_type, post_type)] = avg

    return results


if __name__ == '__main__':
    app = pg.mkQApp()
    pg.dbg()

    all_expts = cached_experiments()
    
    # results = pulse_average_matrix(all_expts, clamp_mode='ic', min_duration=25e-3, pulse_ids=[0, 8])

    # 1. collect a list of all experiments containing connections
    conn_expts = {}
    for c in all_expts.connection_summary():
        conn_expts.setdefault(c['expt'], []).append(c['cells'])

    # 2. for each experiment, get and cache the full set of first-pulse averages
    types = ['sim1', 'tlx3', 'pvalb', 'sst', 'vip']
    indplots = PlotGrid()
    indplots.set_shape(len(types), len(types))
    indplots.show()

    cachefile = "first_pulse_average.pkl"
    cache = pickle.load(open(cachefile, 'rb')) if os.path.isfile(cachefile) else {}

    for expt, conns in sorted(conn_expts.items()):
        if expt.source_id not in cache:
            print "Load:", expt.source_id
            try:
                with expt.data:
                    analyzer = MultiPatchExperimentAnalyzer.get(expt.data)
                    responses = []
                    for pre_cell, post_cell in conns:
                        pre_id, post_id = pre_cell.cell_id-1, post_cell.cell_id-1
                        resp = analyzer.get_evoked_responses(pre_id, post_id, 
                            clamp_mode='ic', 
                            min_duration=25e-3, 
                            pulse_ids=[1]
                        )
                        avg = resp.bsub_mean()
                        cretypes = (pre_cell.cre_type, post_cell.cre_type)

                        if avg is not None:
                            avg.t0 = 0
                            responses.append({'types': cretypes, 'channels': (pre_id, post_id), 'response': avg})
                    cache[expt.source_id] = responses
            except IOError:
                print "Skipped:", expt
                continue

            pkl = pickle.dumps(cache)
            open(cachefile, 'wb').write(pkl)


        for conn in cache[expt.source_id]:
            cretypes = conn['types']
            avg = conn['response']
            if cretypes[0] not in types or cretypes[1] not in types:
                continue
            plt = indplots[types.index(cretypes[0]), types.index(cretypes[1])]
            plt.plot(avg.time_values, avg.data)
            plt.setLabels(left=('%s-%s'%cretypes, 'V'), bottom=('Time', 's'))
            app.processEvents()
            
        
        mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if mem > 12000000:
            print("It is a good day to die.")
            sys.exit(-1)


    # 3. plot matrix of averages by cell type
    plots = PlotGrid()
    plots.set_shape(len(types), len(types))
    plots.show()

    amps = OrderedDict()
    for i, pre_type in enumerate(types):
        for j, post_type in enumerate(types):

            # collect average 
            rg = EvokedResponseGroup(None, None)
            for expt,conns in cache.items():
                for conn in conns:
                    if conn['types'] != (pre_type, post_type):
                        continue
                    if conn['response'].data.std() > 1e-3:
                        print "rejected:", expt, conn
                        continue
                    rg.add(conn['response'], None)

            avg = rg.mean()
            if avg is None:
                continue

            plt = plots[i, j]
            plt.setLabels(left=('%s-%s  n=%d'%(pre_type, post_type, len(rg)), 'V'), bottom=('Time', 's'))
            base = np.median(avg.data[:100])
            plt.plot(avg.time_values, avg.data - base)

            # for first pulse
            fit = fit_psp(avg, yoffset=0, amp_ratio=(0, 'fixed'), mask_stim_artifact=True)
            # for later pulses
            #fit = fit_psp(avg, yoffset=0, mask_stim_artifact=True)
            
            amps[(pre_type, post_type)] = fit.best_values['amp']
            plt.plot(avg.time_values, fit.eval()-base, pen='g')


    plots.setXLink(plots[0,0])
    for i in range(2):
        for j in range(len(types)):
            plots[i,j].setYLink(plots[0,0])
    for i in range(2,len(types)):
        for j in range(len(types)):
            plots[i,j].setYLink(plots[2,0])


    plt = pg.plot()
    typs = amps.keys()
    bar = pg.BarGraphItem(x=np.arange(len(amps)), width=0.6, height=np.array(amps.values()))
    plt.addItem(bar)

    ax = plt.getAxis('bottom')
    ax.setTicks([[(i,"%s-%s"%amps.keys()[i]) for i in range(len(amps))]])


