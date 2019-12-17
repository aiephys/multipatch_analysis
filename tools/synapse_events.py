import pyqtgraph as pg
from neuroanalysis.data import TSeries
from neuroanalysis.fitting import StackedPsp
from aisynphys.database import default_db as db
from aisynphys.data import PulseResponseList
from aisynphys.ui.experiment_browser import ExperimentBrowser


class SynapseEventWindow(pg.QtGui.QSplitter):
    def __init__(self):
        self.loaded_pair = None
        
        pg.QtGui.QSplitter.__init__(self, pg.QtCore.Qt.Horizontal)
        self.ctrl_split = pg.QtGui.QSplitter(pg.QtCore.Qt.Vertical)
        self.addWidget(self.ctrl_split)
        
        self.browser = ExperimentBrowser()
        self.ctrl_split.addWidget(self.browser)
        
        self.scatter_plot = pg.ScatterPlotWidget()
        self.ctrl_split.addWidget(self.scatter_plot.ctrlPanel)
        self.addWidget(self.scatter_plot.plot)

        # set up scatter plot fields
        fields = [
            ('clamp_mode', {'mode': 'enum', 'values': ['ic', 'vc']}),
            ('pulse_number', {'mode': 'range'}),
            ('induction_frequency', {'mode': 'range'}),
            ('recovery_delay', {'mode': 'range'}),
        ]
        fit_keys = ['amp', 'rise_time', 'decay_tau', 'exp_amp', 'yoffset', 'latency', 'nrmse']
        fields = fields + [('fit_'+key, {'mode': 'range'}) for key in fit_keys]
        fields = fields + [('baseline_fit_'+key, {'mode': 'range'}) for key in fit_keys]        
        self.scatter_plot.setFields(fields)
        
        # default filter for IC data
        cm_filter = self.scatter_plot.filter.addNew('clamp_mode')
        cm_filter['vc'] = False
        
        # default color by nrmse 
        nrmse_color = self.scatter_plot.colorMap.addNew('fit_nrmse')

        self.view = pg.GraphicsLayoutWidget()
        self.addWidget(self.view)
        
        self.spike_plot = self.view.addPlot()
        self.data_plot = self.view.addPlot(row=1, col=0)
        self.spike_plot.setXLink(self.data_plot)
        
        self.resize(1600, 1000)
        
        self.browser.itemSelectionChanged.connect(self.browser_item_selected)
        self.scatter_plot.sigScatterPlotClicked.connect(self.scatter_plot_clicked)

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
            self.loaded_pair = pair

            # Load data for scatter plot            
            q = db.query(db.PulseResponse.id.label('prid'), db.PulseResponseFit, db.PatchClampRecording.clamp_mode, db.StimPulse.pulse_number, db.MultiPatchProbe.induction_frequency, db.MultiPatchProbe.recovery_delay)
            q = q.join(db.PulseResponseFit)
            q = q.join(db.StimPulse, db.PulseResponse.stim_pulse)
            q = q.join(db.Recording, db.PulseResponse.recording)
            q = q.join(db.PatchClampRecording)
            q = q.join(db.MultiPatchProbe)
            q = q.filter(db.PulseResponse.pair_id==pair.id)
            df = q.dataframe()
                        
            self.scatter_plot.setData(df.to_records())
            
    def scatter_plot_clicked(self, plt, points):
        self.data_plot.clear()
        self.spike_plot.clear()
        psp = StackedPsp()

        # query raw data for selected points
        ids = [int(pt.data()['pulse_response_id']) for pt in points]
        q = db.query(db.PulseResponse, db.PulseResponse.data, db.PulseResponseFit).join(db.PulseResponseFit).filter(db.PulseResponse.id.in_(ids))
        recs = q.all()
        
        # get pre- and postsynaptic data spike aligned
        prl = PulseResponseList([rec.PulseResponse for rec in recs])
        pre_tsl = prl.pre_tseries(align='spike')
        post_tsl = prl.post_tseries(align='spike')
        
        # plot all
        for i in range(len(prl)):
            pre_ts = pre_tsl[i]
            post_ts = post_tsl[i]
            
            self.spike_plot.plot(pre_ts.time_values, pre_ts.data)
            self.data_plot.plot(post_ts.time_values, post_ts.data)

            # evaluate recorded fit for this response
            fit_par = recs[i].PulseResponseFit
            if fit_par.fit_amp is None:
                continue
            fit = psp.eval(
                x=post_ts.time_values, 
                exp_amp=fit_par.fit_exp_amp,
                exp_tau=fit_par.fit_decay_tau,
                amp=fit_par.fit_amp,
                rise_time=fit_par.fit_rise_time,
                decay_tau=fit_par.fit_decay_tau,
                xoffset=fit_par.fit_latency,
                yoffset=fit_par.fit_yoffset,
                rise_power=2,
            )
            self.data_plot.plot(post_ts.time_values, fit, pen=(0, 255, 0, 100))
            
        self.scatter_plot.setSelectedPoints(points)
        


if __name__ == '__main__':
    import sys
        
    app = pg.mkQApp()
    if sys.flags.interactive == 1:
        pg.dbg()
    
    win = SynapseEventWindow()
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
