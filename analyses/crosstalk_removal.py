"""
Question: can we digitally remove pipette electrical / capacitive crosstalk somehow?

"""
import numpy as np
import pyqtgraph as pg

from neuroanalysis.data import TraceList
from neuroanalysis import filter
from multipatch_analysis.ui.experiment_browser import ExperimentBrowser
from multipatch_analysis.connection_strength import get_amps, get_baseline_amps
from multipatch_analysis import database as db
from multipatch_analysis import data


class StimResponseList(object):
    """A collection of stimulus-response objects that provides some common analysis tools
    """
    def __init__(self, srs):
        self.srs = srs
        self._plot_items = []
        
    def clear_plots(self):
        for item, plt in self._plot_items:
            plt.removeItem(item)
        self._plot_items = []

    def get_tseries(self, series, bsub=True, align='stim', bsub_window=(-3e-3, 0)):
        """Return a TraceList of timeseries, optionally baseline-subtracted and time-aligned.
        
        Parameters
        ----------
        series : str
            "stim", "pre", or "post"
        """
        assert series in ('stim', 'pre', 'post'), "series must be one of 'stim', 'pre', or 'post'"
        tseries = []
        for i,sr in enumerate(self.srs):
            ts = getattr(sr, series + '_tseries')
            if bsub:
                bstart = sr.stim_pulse.onset_time + bsub_window[0]
                bstop = sr.stim_pulse.onset_time + bsub_window[1]
                baseline = np.median(ts.time_slice(bstart, bstop).data)
                ts = ts - baseline
            if align is not None:
                if align == 'stim':
                    t_align = sr.stim_pulse.onset_time 
                elif align == 'pre':
                    t_align = sr.stim_pulse.spikes[0].max_dvdt_time
                elif align == 'post':
                    raise NotImplementedError()
                else:
                    raise ValueError("invalid time alignment mode %r" % align)
                ts = ts.copy(t0=ts.t0-t_align)
            tseries.append(ts)
        return TraceList(tseries)

    def plot_stimulus(self, plt, **kwds):
        stim_ts = self.get_tseries('stim', **kwds)
        self._plot_ts(stim_ts, plt)
    
    def plot_presynaptic(self, plt, **kwds):
        pre_ts = self.get_tseries('pre', **kwds)
        self._plot_ts(pre_ts, plt)

    def plot_postsynaptic(self, plt, **kwds):
        post_ts = self.get_tseries('post', **kwds)
        self._plot_ts(post_ts, plt)
        
    def _plot_ts(self, ts_list, plt):        
        for ts in ts_list:
            item = plt.plot(ts.time_values, ts.data, pen=(255, 255, 255, 80), antialias=True)
            self._plot_items.append((item, plt))
        avg = ts_list.mean()
        item = plt.plot(avg.time_values, avg.data, pen={'color':'g', 'width':2}, shadowPen={'color':'k', 'width':3}, antialias=True)
        self._plot_items.append((item, plt))

    def __iter__(self):
        for sr in self.srs:
            yield sr

    def __getitem__(self, index):
        item = self.srs.__getitem__(index)
        if isinstance(item, list):
            return StimResponseList(item)
        else:
            return item


class CrosstalkAnalyzer(object):
    """User interface for exploring crosstalk removal methods
    """
    def __init__(self, expt_browser=None):
        if expt_browser is None:
            expt_browser = ExperimentBrowser()
        self.expt_browser = expt_browser
        expt_browser.itemSelectionChanged.connect(self.expt_selection_changed)
        self.session = db.Session()        
        self.response_list = None
        self.corrected_response_list = None

        self.params = pg.parametertree.Parameter(name="crosstalk removal", type="group", children=[
            dict(name='limit', type='int', value=20),
            dict(name='correction', type='group', children=[
                dict(name='remove stim', type='bool', value=True),
                dict(name='stim scale', type='float', value=30e3, dec=True, step=0.5),
                dict(name='stim lowpass', type='float', value=7e3, dec=True, step=0.5, suffix='Hz', siPrefix=True),
                dict(name='remove spike', type='bool', value=False),
                dict(name='spike scale', type='float', value=0.5, dec=True, step=0.5),
                dict(name='spike lowpass', type='float', value=10e3, dec=True, step=0.5, suffix='Hz', siPrefix=True),
                dict(name='remove spike dv/dt', type='bool', value=True),
                dict(name='spike dv/dt scale', type='float', value=0.03, dec=True, step=0.5),
                dict(name='spike dv/dt lowpass', type='float', value=550., dec=True, step=0.5, suffix='Hz', siPrefix=True),
            ]),
            dict(name='display', type='group', children=[
                dict(name='plot spike aligned', type='bool', value=False),
                dict(name='plot dv/dt', type='bool', value=False),
                dict(name='plot lowpass', type='float', value=10e3, dec=True, step=0.5, suffix='Hz', siPrefix=True),
            ])
        ])
        self.ptree = pg.parametertree.ParameterTree()
        self.ptree.setParameters(self.params)
        self.params.child('correction').sigTreeStateChanged.connect(self.update_analysis)
        self.params.child('display').sigTreeStateChanged.connect(self.update_display)

        self.pw = pg.GraphicsLayoutWidget()
        
        self.plt1 = self.pw.addPlot(0, 0)
        self.plt2 = self.pw.addPlot(1, 0)
        self.plt3 = self.pw.addPlot(2, 0)
        self.plt2.setXLink(self.plt1)
        self.plt3.setXLink(self.plt1)

    def show(self):
        self.win = pg.QtGui.QSplitter(pg.QtCore.Qt.Horizontal)
        self.vsplit = pg.QtGui.QSplitter(pg.QtCore.Qt.Vertical)
        self.win.addWidget(self.vsplit)
        self.vsplit.addWidget(self.expt_browser)
        self.vsplit.addWidget(self.ptree)
        self.win.addWidget(self.pw)
        self.win.show()
        return self.win
        
    def expt_selection_changed(self):
        sel = self.expt_browser.selectedItems()
        if len(sel) == 0:
            return
        sel = sel[0]
        if not hasattr(sel, 'pair'):
            return
        self.load_pair(sel.pair)

    def select_pair(self, pair_id):
        self.expt_browser.select_pair(pair_id)

    def load_pair(self, pair):
        with pg.BusyCursor():
            # preload pulse response data
            db.query(db.PulseResponse.data).filter(db.PulseResponse.pair_id==pair.id).all()
            
            self.response_list = StimResponseList(pair.pulse_responses)
            self.update_analysis()

    def update_analysis(self):
        if self.corrected_response_list is not None:
            self.corrected_response_list.clear_plots()
        
        limit = self.params['limit']
        selected_responses = self.response_list[:limit] if limit > 0 else self.response_list
        self.corrected_response_list = StimResponseList([self.process_response(sr) for sr in selected_responses])
        self.update_display()

    def update_display(self):
        kwds = {
            'align': 'pre' if self.params['display', 'plot spike aligned'] else 'stim',
        } 
        self.corrected_response_list.clear_plots()
        self.corrected_response_list.plot_stimulus(self.plt1, **kwds)
        self.corrected_response_list.plot_presynaptic(self.plt2, **kwds)
        self.corrected_response_list.plot_postsynaptic(self.plt3, **kwds)

    def process_response(self, sr):
        stim = sr.stim_tseries
        pre = sr.pre_tseries
        post = sr.post_tseries.copy()

        if self.params['correction', 'remove stim']:
            pulse_start = sr.stim_pulse.onset_time
            pulse_stop = pulse_start + sr.stim_pulse.duration
            stim_correction = post.copy(data=np.zeros(len(post)))
            stim_correction.time_slice(pulse_start, pulse_stop).data[:] = self.params['correction', 'stim scale'] * sr.stim_pulse.amplitude
            stim_correction = filter.bessel_filter(stim_correction, self.params['correction', 'stim lowpass'], bidir=False)
            post = post - stim_correction.data

        if self.params['correction', 'remove spike']:
            spike = pre.copy(data=pre.data * self.params['correction', 'spike scale'])
            spike_correction = filter.bessel_filter(spike, self.params['correction', 'spike lowpass'], bidir=False)
            start = max(spike.t0, post.t0)
            stop = min(spike.t_end, post.t_end)

            post = post.time_slice(start, stop) - spike_correction.time_slice(start, stop).data

        if self.params['correction', 'remove spike dv/dt']:
            spike_diff = np.diff(pre.data) * self.params['correction', 'spike dv/dt scale']
            spike = pre.copy(data=spike_diff)
            spike_correction = filter.bessel_filter(spike, self.params['correction', 'spike dv/dt lowpass'], bidir=False)
            start = max(spike.t0, post.t0)
            stop = min(spike.t_end, post.t_end)

            post = post.time_slice(start, stop) - spike_correction.time_slice(start, stop).data

        return data.PulseResponse(stim_pulse=sr.stim_pulse, post_tseries=post)

        

if __name__ == '__main__':
    import sys
    app = pg.mkQApp()
    pg.dbg()
    
    cta = CrosstalkAnalyzer()

    if len(sys.argv) > 1:
        pair_id = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
        cta.select_pair(pair_id)

    win = cta.show()
    win.resize(1400, 1000)
    
    if sys.flags.interactive == 0:
        app.exec_()
