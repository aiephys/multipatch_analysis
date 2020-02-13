import pandas as pd
import pyqtgraph as pg
from neuroanalysis.data import TSeries
from neuroanalysis.fitting import StackedPsp, Psp
from aisynphys.database import default_db as db
from aisynphys.data import PulseResponseList
from aisynphys.ui.experiment_browser import ExperimentBrowser
from aisynphys.pulse_response_strength import deconv_filter
from aisynphys import config


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
            ('ex_qc_pass', {'mode': 'enum', 'values': [True, False]}),
            ('in_qc_pass', {'mode': 'enum', 'values': [True, False]}),
            ('pulse_number', {'mode': 'enum', 'values': list(range(1,13))}),
            ('induction_frequency', {'mode': 'range'}),
            ('recovery_delay', {'mode': 'range'}),
            ('pos_amp', {'mode': 'range'}),
            ('neg_amp', {'mode': 'range'}),
            ('pos_dec_amp', {'mode': 'range'}),
            ('neg_dec_amp', {'mode': 'range'}),
        ]
        fit_keys = ['amp', 'rise_time', 'decay_tau', 'exp_amp', 'yoffset', 'latency', 'nrmse']
        fields = fields + [('fit_'+key, {'mode': 'range'}) for key in fit_keys]
        fields = fields + [('baseline_fit_'+key, {'mode': 'range'}) for key in fit_keys]
        fields = fields + [('dec_fit_'+key, {'mode': 'range'}) for key in fit_keys]
        fields = fields + [('baseline_dec_fit_'+key, {'mode': 'range'}) for key in fit_keys]
        fields = fields + [
            ('dec_fit_reconv_amp', {'mode': 'range'}),
            ('baseline_dec_fit_reconv_amp', {'mode': 'range'}),
        ]
        self.scatter_plot.setFields(fields)
        
        # default filter for IC data
        cm_filter = self.scatter_plot.filter.addNew('clamp_mode')
        cm_filter['vc'] = False
        
        # default color by nrmse 
        nrmse_color = self.scatter_plot.colorMap.addNew('fit_nrmse')
        qc_color = self.scatter_plot.colorMap.addNew('ex_qc_pass')
        qc_color['Values', 'True'] = (255, 255, 255)
        qc_color['Values', 'False'] = (0, 0, 0)
        qc_color['Operation'] = 'Multiply'
        

        self.view = pg.GraphicsLayoutWidget()
        self.addWidget(self.view)
        
        self.spike_plots = [self.view.addPlot(row=0, col=0), self.view.addPlot(row=0, col=1)]
        self.data_plots = [self.view.addPlot(row=1, col=0), self.view.addPlot(row=1, col=1)]
        self.dec_plots = [self.view.addPlot(row=2, col=0), self.view.addPlot(row=2, col=1)]
        for col in (0, 1):
            self.spike_plots[col].setXLink(self.data_plots[0])
            self.data_plots[col].setXLink(self.data_plots[0])
            self.dec_plots[col].setXLink(self.data_plots[0])

        self.spike_plots[1].setYLink(self.spike_plots[0])
        self.data_plots[1].setYLink(self.data_plots[0])
        self.dec_plots[1].setYLink(self.dec_plots[0])
        
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
            q = db.query(
                db.PulseResponse.id.label('prid'), 
                db.PulseResponse.ex_qc_pass,
                db.PulseResponse.in_qc_pass,
                db.PulseResponseFit, 
                
                db.PulseResponseStrength.pos_amp,
                db.PulseResponseStrength.neg_amp,
                db.PulseResponseStrength.pos_dec_amp,
                db.PulseResponseStrength.neg_dec_amp,
                db.PulseResponseStrength.pos_dec_latency,
                db.PulseResponseStrength.neg_dec_latency,
                db.PulseResponseStrength.crosstalk,

                # db.BaselineResponseStrength,
                db.PatchClampRecording.clamp_mode,
                db.StimPulse.pulse_number, 
                db.MultiPatchProbe.induction_frequency, 
                db.MultiPatchProbe.recovery_delay
            )
            q = q.join(db.PulseResponseFit)
            q = q.join(db.PulseResponseStrength)
            # q = q.join(db.BaselineResponseStrength)
            q = q.join(db.StimPulse, db.PulseResponse.stim_pulse)
            q = q.join(db.Recording, db.PulseResponse.recording)
            q = q.join(db.PatchClampRecording)
            q = q.join(db.MultiPatchProbe)
            q = q.filter(db.PulseResponse.pair_id==pair.id)
            
            df = q.dataframe()
            for col in df.columns:
                if 'fit_' in col:
                    df[col] = pd.to_numeric(df[col])
            
            self.scatter_plot.setData(df.to_records())
            
    def scatter_plot_clicked(self, plt, points):
        for plt in self.data_plots + self.dec_plots + self.spike_plots:
            plt.clear()

        # query raw data for selected points
        ids = [int(pt.data()['pulse_response_id']) for pt in points]
        q = db.query(db.PulseResponse, db.PulseResponse.data, db.PulseResponseFit).join(db.PulseResponseFit).filter(db.PulseResponse.id.in_(ids))
        recs = q.all()
        
        for rec in recs:
            self._plot_pulse_response(rec)
            
        self.scatter_plot.setSelectedPoints(points)
        

    def _plot_pulse_response(self, rec):
        pr = rec.PulseResponse
        
        pre_ts = pr.get_tseries('pre', align_to='spike')
        
        # If there is no presynaptic spike time, plot spike in red and bail out
        if pre_ts is None:
            pre_ts = pr.get_tseries('pre', align_to='pulse')
            self.spike_plots[0].plot(pre_ts.time_values, pre_ts.data, pen=(255, 0, 0, 100))
            return
        
        post_ts = pr.get_tseries('post', align_to='spike')
        
        self.spike_plots[0].plot(pre_ts.time_values, pre_ts.data)
        self.data_plots[0].plot(post_ts.time_values, post_ts.data)

        # evaluate recorded fit for this response
        fit_par = rec.PulseResponseFit
        
        # If there is no fit, bail out here
        if fit_par.fit_amp is None:
            return
        
        spsp = StackedPsp()
        fit = spsp.eval(
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
        self.data_plots[0].plot(post_ts.time_values, fit, pen=(0, 255, 0, 100))

        # plot with reconvolved amplitude
        fit = spsp.eval(
            x=post_ts.time_values, 
            exp_amp=fit_par.fit_exp_amp,
            exp_tau=fit_par.fit_decay_tau,
            amp=fit_par.dec_fit_reconv_amp,
            rise_time=fit_par.fit_rise_time,
            decay_tau=fit_par.fit_decay_tau,
            xoffset=fit_par.fit_latency,
            yoffset=fit_par.fit_yoffset,
            rise_power=2,
        )
        self.data_plots[0].plot(post_ts.time_values, fit, pen=(200, 255, 0, 100))

        # plot deconvolution
        clamp_mode = pr.recording.patch_clamp_recording.clamp_mode
        if clamp_mode == 'ic':
            decay_tau = self.loaded_pair.synapse.psp_decay_tau
            lowpass = 2000
        else:
            decay_tau = self.loaded_pair.synapse.psc_decay_tau
            lowpass = 6000
            
        dec = deconv_filter(post_ts, None, tau=decay_tau, lowpass=lowpass, remove_artifacts=False, bsub=True)
        self.dec_plots[0].plot(dec.time_values, dec.data)

        # plot deconvolution fit
        psp = Psp()
        fit = psp.eval(
            x=dec.time_values,
            exp_tau=fit_par.dec_fit_decay_tau,
            amp=fit_par.dec_fit_amp,
            rise_time=fit_par.dec_fit_rise_time,
            decay_tau=fit_par.dec_fit_decay_tau,
            xoffset=fit_par.dec_fit_latency,
            yoffset=fit_par.dec_fit_yoffset,
            rise_power=1,
        )
        self.dec_plots[0].plot(dec.time_values, fit, pen=(0, 255, 0, 100))


if __name__ == '__main__':
    import sys, argparse
    
    parser = argparse.ArgumentParser(parents=[config.parser])
    parser.add_argument('experiment_id', type=str, nargs='?')
    parser.add_argument('pre_cell_id', type=str, nargs='?')
    parser.add_argument('post_cell_id', type=str, nargs='?')
    args = parser.parse_args()
        
    app = pg.mkQApp()
    if sys.flags.interactive == 1:
        pg.dbg()
    
    win = SynapseEventWindow()
    win.show()
    
    if args.post_cell_id is not None:
        expt = db.experiment_from_ext_id(args.experiment_id)
        win.browser.populate([expt], synapses=True)
        pair = expt.pairs[args.pre_cell_id, args.post_cell_id]
        win.browser.select_pair(pair.id)
    else:
        win.browser.populate(synapses=True)

    if sys.flags.interactive == 0:
        app.exec_()
