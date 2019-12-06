import pyqtgraph as pg
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import TSeries
from neuroanalysis.fitting import StackedPsp
from aisynphys.database import default_db as db
from aisynphys.ui.experiment_browser import ExperimentBrowser
from aisynphys.dynamics import pulse_response_query, sorted_pulse_responses


class DynamicsWindow(pg.QtGui.QSplitter):
    def __init__(self):
        self.loaded_pair = None
        
        pg.QtGui.QSplitter.__init__(self, pg.QtCore.Qt.Horizontal)
        self.ctrl_split = pg.QtGui.QSplitter(pg.QtCore.Qt.Vertical)
        self.addWidget(self.ctrl_split)
        
        self.browser = ExperimentBrowser()
        self.ctrl_split.addWidget(self.browser)
        
        self.ptree = pg.parametertree.ParameterTree()
        self.ctrl_split.addWidget(self.ptree)
        
        self.params = pg.parametertree.Parameter.create(name='params', type='group', children=[
            {'name': 'show spikes', 'type': 'bool', 'value': False},
            {'name': 'stimulus filter', 'type': 'group'},
        ])
        self.ptree.setParameters(self.params)
        
        self.view = pg.GraphicsLayoutWidget()
        self.addWidget(self.view)
        
        self.plots = []
        
        self.browser.itemSelectionChanged.connect(self.browser_item_selected)
        self.params.sigTreeStateChanged.connect(self.plot_all)

    def clear(self):
        for plt in self.plots:
            self.view.removeItem(plt)
        self.plots = []
    
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
        if pair is not self.loaded_pair:
            print("Loading:", pair)
            q = pulse_response_query(pair, data=True)
            self.sorted_recs = sorted_pulse_responses(q.all())
            self.stim_keys = sorted(list(self.sorted_recs.keys()))
            self.update_params()
            self.loaded_pair = pair
        
        self.plot_all()
        
    def update_params(self):
        with pg.SignalBlock(self.params.sigTreeStateChanged, self.plot_all):
            stim_param = self.params.child('stimulus filter')
            for ch in stim_param.children():
                stim_param.removeChild(ch)
            
            for k in self.stim_keys:
                param = pg.parametertree.Parameter.create(name=str(k), type="bool", value="True")
                stim_param.addChild(param)
        
    def plot_all(self):
        self.clear()

        for i,stim_key in enumerate(self.stim_keys):
            if self.params['stimulus filter', str(stim_key)] is False:
                continue
            
            plt = DynamicsPlot()
            self.plots.append(plt)
            self.view.addItem(plt)
            self.view.nextRow()
            plt.set_title("%s  %0.0f Hz  %0.2f s" % stim_key)
            prs = self.sorted_recs[stim_key]
            plt.set_data(prs)
        

class DynamicsPlot(pg.GraphicsLayout):
    def __init__(self):
        pg.GraphicsLayout.__init__(self)
        self.show_spikes = False
        
        self.label = pg.TextItem()
        self.label.setParentItem(self)
        self.label.setPos(200, 0)
        self.spike_plot = self.addPlot()
        self.spike_plot.setVisible(False)
        self.nextRow()
        self.data_plot = self.addPlot()
        
    def set_title(self, title):
        self.label.setText(title)
        
    def set_data(self, data):
        psp = StackedPsp()
        
        for recording in data:
            pulses = sorted(list(data[recording].keys()))
            for pulse_n in pulses:
                rec = data[recording][pulse_n]
                # spike-align pulse + offset for pulse number
                spike_t = rec.StimPulse.first_spike_time
                if spike_t is None:
                    spike_t = rec.StimPulse.onset_time + 1e-3
                    
                qc_pass = rec.PulseResponse.in_qc_pass if rec.Synapse.synapse_type == 'in' else rec.PulseResponse.ex_qc_pass
                pen = (255, 255, 255, 100) if qc_pass else (100, 0, 0, 100)
                
                t0 = rec.PulseResponse.data_start_time - spike_t
                ts = TSeries(data=rec.data, t0=t0, sample_rate=db.default_sample_rate)
                c = self.data_plot.plot(ts.time_values, ts.data, pen=pen)
                
                # arrange plots nicely
                shift = (pulse_n * 35e-3 + (30e-3 if pulse_n > 8 else 0), 0)
                c.setPos(*shift)
                
                if not qc_pass:
                    c.setZValue(-10)
                    continue
                    
                # evaluate recorded fit for this response
                fit_par = rec.PulseResponseFit
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
                c = self.data_plot.plot(ts.time_values, fit, pen=(0, 255, 0, 100))
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
        
