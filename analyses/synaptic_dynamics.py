import pyqtgraph as pg
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import TSeries
from neuroanalysis.fitting import StackedPsp
from aisynphys.database import default_db as db
from aisynphys.ui.experiment_browser import ExperimentBrowser
from aisynphys.dynamics import pulse_response_query, sorted_pulse_responses


class DynamicsWindow(pg.QtGui.QSplitter):
    def __init__(self):
        pg.QtGui.QSplitter.__init__(self, pg.QtCore.Qt.Horizontal)
        
        self.browser = ExperimentBrowser()
        self.addWidget(self.browser)
        
        self.plots = PlotGrid()
        self.addWidget(self.plots)
        
        self.browser.itemSelectionChanged.connect(self.browser_item_selected)

    def browser_item_selected(self):
        with pg.BusyCursor():
            selected = self.browser.selectedItems()
            if len(selected) != 1:
                return
            item = selected[0]
            if not hasattr(item, 'pair'):
                return
            pair = item.pair

            self.load_pair(pair)
        
    def load_pair(self, pair):
        print("Loading:", pair)
        
        q = pulse_response_query(pair, data=True)
        self.sorted_recs = sorted_pulse_responses(q.all())
        
        self.plot_all()
        
    def plot_all(self):
        self.plots.clear()
        self.plots.set_shape(len(self.sorted_recs), 1)
        psp = StackedPsp()
        
        stim_keys = sorted(list(self.sorted_recs.keys()))
        for i,stim_key in enumerate(stim_keys):
            prs = self.sorted_recs[stim_key]
            plt = self.plots[i,0]
            plt.setTitle("%s  %0.0f Hz  %0.2f s" % stim_key)
            
            
            for recording in prs:
                pulses = sorted(list(prs[recording].keys()))
                for pulse_n in pulses:
                    rec = prs[recording][pulse_n]
                    # spike-align pulse + offset for pulse number
                    spike_t = rec.stim_pulse.first_spike_time
                    if spike_t is None:
                        spike_t = rec.stim_pulse.onset_time + 1e-3
                        
                    qc_pass = rec.pulse_response.in_qc_pass if rec.synapse.synapse_type == 'in' else rec.pulse_response.ex_qc_pass
                    pen = (255, 255, 255, 100) if qc_pass else (100, 0, 0, 100)
                    
                    t0 = rec.pulse_response.data_start_time - spike_t
                    ts = TSeries(data=rec.data, t0=t0, sample_rate=db.default_sample_rate)
                    c = plt.plot(ts.time_values, ts.data, pen=pen)
                    
                    # arrange plots nicely
                    shift = (pulse_n * 35e-3 + (30e-3 if pulse_n > 8 else 0), 0)
                    c.setPos(*shift)
                    
                    if not qc_pass:
                        c.setZValue(-10)
                        continue
                        
                    # evaluate recorded fit for this response
                    fit_par = rec.pulse_response_fit
                    if fit_par.fit_amp is None:
                        continue
                    fit = psp.eval(
                        x=ts.time_values, 
                        exp_amp=fit_par.fit_exp_amp,
                        exp_tau=fit_par.fit_decay_tau,
                        amp=fit_par.fit_amp,
                        rise_time=fit_par.fit_rise_time,
                        decay_tau=fit_par.fit_decay_tau,
                        xoffset=fit_par.fit_latency,
                        yoffset=fit_par.fit_yoffset,
                        rise_power=2,
                    )
                    c = plt.plot(ts.time_values, fit, pen=(0, 255, 0, 100))
                    c.setZValue(10)
                    c.setPos(*shift)
        


if __name__ == '__main__':
    import sys
        
    app = pg.mkQApp()
    pg.dbg()
    
    win = DynamicsWindow()
    win.show()
    
    if len(sys.argv) > 1:
        expt_id, pre_cell, post_cell = sys.argv[1:]
        expt = db.experiment_from_ext_id(expt_id)
        win.browser.populate([expt])
        pair = expt.pairs[pre_cell, post_cell]
        win.browser.select_pair(pair.id)
    else:
        win.browser.populate()

    if sys.flags.interactive == 0:
        app.exec_()
        
