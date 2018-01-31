import sys
import pyqtgraph as pg
from multipatch_analysis.experiment_list import ExperimentList
from multipatch_analysis.synaptic_dynamics import DynamicsAnalyzer


if __name__ == '__main__':
    app = pg.mkQApp()
    pg.dbg()
    
    expt_ind = sys.argv[1]
    all_expts = ExperimentList(cache='expts_cache.pkl')
    expt = all_expts[expt_ind]

    pre_cell = int(sys.argv[2])
    post_cell = int(sys.argv[3])

    method = 'deconv' if '--deconv' in sys.argv else 'fit'

    analyzer = DynamicsAnalyzer(expt, pre_cell, post_cell, method=method)
    if len(analyzer.pulse_responses) == 0:
        raise Exception("No suitable data found for cell %d -> cell %d in expt %s" % (pre_cell, post_cell, expt_ind))
           
    # Plot all individual and averaged train responses for all sets of stimulus parameters
    train_plots = analyzer.plot_train_responses()

    if '--no-fit' in sys.argv:
        sys.exit(0)  # user requested no fitting; bail out early

    if '--deconv' in sys.argv:
        # get deconvolved response trains
        analyzer.plot_deconvolved_trains(train_plots)
        analyzer.measure_train_amps_from_deconv(plot_grid=train_plots)
        
    else:
        # Estimate PSP amplitude
        amp_est, amp_sign, avg_amp, amp_plot, n_sweeps = analyzer.estimate_amplitude(plot=True)
        app.processEvents()
        
        # Estimate PSP kinetics
        rise_time, decay_tau, latency, kin_plot = analyzer.estimate_kinetics(plot=True)
        app.processEvents()
        
        # Fit trains to multi-event models and plot the results
        analyzer.plot_train_fits(train_plots)

    if '--no-model' in sys.argv:
        sys.exit(0)  # user requested no model; bail out early
        
    # update GUI before doing model fit
    app.processEvents()
    
    # Fit release model to dynamics
    model, fit = analyzer.fit_release_model(dynamics=['Dep', 'Fac', 'UR'])
    
    # plot nost recently-generated fit results
    rel_plots = analyzer.plot_model_results()
