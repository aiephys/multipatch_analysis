# coding: utf8
"""
Big question: what's the best way to measure synaptic strength / connectivity?

1. Whatever method we use to measure connectivity, we also need to characterize the detection limit (per synapse)
2. Any method we choose should be run on both pulse response and background data, and the distributions of
   these results must be compared to make the connectivity call

"""
from __future__ import print_function, division

import argparse, time, sys, os, pickle, io
import numpy as np
import scipy.stats
import pandas
from datetime import datetime

from sqlalchemy.orm import aliased
import sklearn.svm, sklearn.preprocessing, sklearn.ensemble

import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore

from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.data import Trace, TraceList
from neuroanalysis import filter
from neuroanalysis.event_detection import exp_deconvolve
from neuroanalysis.baseline import float_mode
from neuroanalysis.fitting import Psp

from multipatch_analysis.database import database as db
from multipatch_analysis.ui.multipatch_nwb_viewer import MultipatchNwbViewer
from multipatch_analysis.pulse_response_strength import (
    PulseResponseStrength, BaselineResponseStrength, response_query,
    baseline_query, analyze_response_strength, pulse_response_strength_tables,
)
from multipatch_analysis.connection_strength import ConnectionStrength, get_amps, get_baseline_amps
import multipatch_analysis.morphology  # just to initialize ORM
from multipatch_analysis import constants


class ExperimentBrowser(pg.TreeWidget):
    def __init__(self):
        pg.TreeWidget.__init__(self)
        self.setColumnCount(7)
        self.setHeaderLabels(['date', 'timestamp', 'rig', 'organism', 'region', 'genotype', 'acsf'])
        self.setDragDropMode(self.NoDragDrop)
        self.populate()
        self._last_expanded = None
        
    def populate(self):
        self.items_by_pair_id = {}
        
        self.session = db.Session()
        db_expts = db.list_experiments(session=self.session)
        db_expts.sort(key=lambda e: e.acq_timestamp)
        for expt in db_expts:
            date = expt.acq_timestamp
            date_str = datetime.fromtimestamp(date).strftime('%Y-%m-%d')
            slice = expt.slice
            expt_item = pg.TreeWidgetItem(map(str, [date_str, '%0.3f'%expt.acq_timestamp, expt.rig_name, slice.species, expt.target_region, slice.genotype, expt.acsf]))
            expt_item.expt = expt
            self.addTopLevelItem(expt_item)

            for pair in expt.pairs:
                if pair.n_ex_test_spikes == 0 and pair.n_in_test_spikes == 0:
                    continue
                cells = '%d => %d' % (pair.pre_cell.ext_id, pair.post_cell.ext_id)
                conn = {True:"syn", False:"-", None:"?"}[pair.synapse]
                types = 'L%s %s => L%s %s' % (pair.pre_cell.target_layer or "?", pair.pre_cell.cre_type, pair.post_cell.target_layer or "?", pair.post_cell.cre_type)
                pair_item = pg.TreeWidgetItem([cells, conn, types])
                expt_item.addChild(pair_item)
                pair_item.pair = pair
                pair_item.expt = expt
                self.items_by_pair_id[pair.id] = pair_item
                
    def select_pair(self, pair_id):
        """Select a specific pair from the list
        """
        if self._last_expanded is not None:
            self._last_expanded.setExpanded(False)
        item = self.items_by_pair_id[pair_id]
        self.clearSelection()
        item.setSelected(True)
        parent = item.parent()
        if not parent.isExpanded():
            parent.setExpanded(True)
            self._last_expanded = parent
        else:
            self._last_expanded = None
        self.scrollToItem(item)





class ResponseStrengthPlots(pg.dockarea.DockArea):
    def __init__(self, session):
        pg.dockarea.DockArea.__init__(self)
        self.session = session

        self.analyses = [('neg', 'ic'), ('pos', 'ic'), ('neg', 'vc'), ('pos', 'vc')]
        self.analyzers = []
        self.analyzer_docks = []
        for col, analysis in enumerate(self.analyses):
            analyzer = ResponseStrengthAnalyzer(analysis, session)
            self.analyzers.append(analyzer)
            dock = pg.dockarea.Dock(analyzer.title, widget=analyzer.widget)
            self.analyzer_docks.append(dock)
            self.addDock(dock, 'right')
                
    def load_conn(self, pair):
        with pg.BusyCursor():
            for analyzer in self.analyzers:
                analyzer.load_conn(pair)


class ResponseStrengthAnalyzer(object):
    def __init__(self, analysis, session):
        self.analysis = analysis  # ('pos'|'neg', 'ic'|'vc')
        self.title = ' '.join(analysis)
        self.session = session
        self._amp_recs = None
        self._base_recs = None

        self.widget = QtGui.QWidget()
        self.layout = QtGui.QGridLayout()
        self.widget.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)

        self.gl = pg.GraphicsLayoutWidget()
        self.layout.addWidget(self.gl, 0, 0)

        # histogram plots
        self.hist_plot = pg.PlotItem(title=self.title)
        self.gl.addItem(self.hist_plot, row=0, col=0)
        self.hist_plot.hideAxis('bottom')

        # event scatter plots
        self.scatter_plot = pg.PlotItem()
        self.gl.addItem(self.scatter_plot, row=1, col=0)
        self.scatter_plot.getAxis('left').setTicks([[(1.4, 'fg'), (0.4, 'bg')], []])
        self.scatter_plot.setFixedHeight(200)
        self.hist_plot.setXLink(self.scatter_plot)

        self.bg_scatter = pg.ScatterPlotItem(pen=None, symbol='o', symbolPen=None)
        self.fg_scatter = pg.ScatterPlotItem(pen=None, symbol='o', symbolPen=None)
        self.scatter_plot.addItem(self.bg_scatter)
        self.scatter_plot.addItem(self.fg_scatter)
        self.fg_scatter.sigClicked.connect(self.fg_scatter_clicked)
        self.bg_scatter.sigClicked.connect(self.bg_scatter_clicked)

        # event selector
        self.event_region = pg.LinearRegionItem()
        self.scatter_plot.addItem(self.event_region)
        self.event_region.setZValue(-10)
        self.event_region.sigRegionChangeFinished.connect(self.region_changed)

        # trace plots
        self.fg_trace_plot = pg.PlotItem()
        self.gl.addItem(self.fg_trace_plot, row=2, col=0)
        self.fg_trace_plot.hideAxis('bottom')
        self.fg_trace_plot.addLine(x=10e-3, pen=(0, 0, 100), movable=False)
        self.bg_trace_plot = pg.PlotItem(labels={'bottom': ('time', 's')})
        self.gl.addItem(self.bg_trace_plot, row=3, col=0)
        self.bg_trace_plot.setXLink(self.fg_trace_plot)
        self.bg_trace_plot.setYLink(self.fg_trace_plot)
        self.bg_trace_plot.addLine(x=10e-3, pen=(0, 0, 100), movable=False)
        self.fg_trace_plot.setXRange(0, 20e-3)

        self.ctrl = QtGui.QWidget()
        self.layout.addWidget(self.ctrl, 1, 0)
        self.ctrl_layout = QtGui.QGridLayout()
        self.ctrl.setLayout(self.ctrl_layout)
        self.ctrl_layout.setContentsMargins(0, 0, 0, 0)

        self.field_combo = QtGui.QComboBox()
        for field in ['dec_amp', 'amp', 'dec_latency', 'crosstalk']:
            self.field_combo.addItem(field)
        self.ctrl_layout.addWidget(self.field_combo, 0, 0)
        self.field_combo.currentIndexChanged.connect(self.update_scatter_plots)

        self.deconv_check = QtGui.QCheckBox('deconvolve')
        self.deconv_check.setChecked(True)
        self.ctrl_layout.addWidget(self.deconv_check, 1, 0)
        self.deconv_check.toggled.connect(self.replot_all)

        self.bsub_check = QtGui.QCheckBox('bsub')
        self.bsub_check.setChecked(True)
        self.ctrl_layout.addWidget(self.bsub_check, 1, 1)
        self.bsub_check.toggled.connect(self.replot_all)

        self.lpf_check = QtGui.QCheckBox('lpf')
        self.lpf_check.setChecked(True)
        self.ctrl_layout.addWidget(self.lpf_check, 1, 2)
        self.lpf_check.toggled.connect(self.replot_all)

        self.ar_check = QtGui.QCheckBox('crosstalk')
        self.ar_check.setChecked(True)
        self.ctrl_layout.addWidget(self.ar_check, 1, 3)
        self.ar_check.toggled.connect(self.replot_all)

        self.align_check = QtGui.QCheckBox('align')
        self.align_check.setChecked(True)
        self.ctrl_layout.addWidget(self.align_check, 1, 4)
        self.align_check.toggled.connect(self.replot_all)

        self.pulse_ctrl = QtGui.QWidget()
        self.ctrl_layout.addWidget(self.pulse_ctrl, 2, 0, 1, 5)
        self.pulse_layout = QtGui.QHBoxLayout()
        self.pulse_ctrl.setLayout(self.pulse_layout)
        self.pulse_layout.setContentsMargins(0, 0, 0, 0)
        self.pulse_layout.setSpacing(0)

        self.color_by_pulse_check = QtGui.QCheckBox('color pulse n')
        self.pulse_layout.addWidget(self.color_by_pulse_check)
        self.color_by_pulse_check.toggled.connect(self.update_scatter_plots)

        self.pulse_checks = []
        for i in range(12):
            c = QtGui.QCheckBox()
            c.setChecked(True)
            self.pulse_checks.append(c)
            self.pulse_layout.addWidget(c)
            c.setMaximumWidth(20)
            c.toggled.connect(self.update_scatter_plots)

        self.selected_fg_traces = []
        self.selected_bg_traces = []
        self.clicked_fg_traces = []
        self.clicked_bg_traces = []
        self._selected_fg_ids = []
        self._selected_bg_ids = []
        self._clicked_fg_ids = []
        self._clicked_bg_ids = []
    
    def fg_scatter_clicked(self, sp, points):
        """Point(s) were clicked; plot their source traces in a different color.
        """
        ids = [p.data()['id'] for p in points]
        self._clicked_fg_ids = ids
        self.plot_prd_ids(ids, 'fg', trace_list=self.clicked_fg_traces, pen='y')

        global selected_response
        selected_response = self.session.query(PulseResponseStrength).filter(PulseResponseStrength.id==ids[0]).first()
        prs_qc()

    def bg_scatter_clicked(self, sp, points):
        """Point(s) were clicked; plot their source traces in a different color.
        """
        ids = [p.data()['id'] for p in points]
        self._clicked_bg_ids = ids
        self.plot_prd_ids(ids, 'bg', trace_list=self.clicked_bg_traces, pen='y')

        global selected_response
        selected_response = self.session.query(BaselineResponseStrength).filter(BaselineResponseStrength.id==ids[0]).first()
        prs_qc()

    def load_conn(self, pair):
        self.pair = pair
        self._amp_recs = get_amps(self.session, self.pair, clamp_mode=self.analysis[1])
        self._base_recs = get_baseline_amps(self.session, self.pair, amps=self._amp_recs, clamp_mode=self.analysis[1])
        self.update_scatter_plots()

    def update_scatter_plots(self): 
        amp_recs = self._amp_recs
        base_recs = self._base_recs
        if amp_recs is None:
            return

        # select fg/bg data
        fg_data = amp_recs
        bg_data = base_recs[:len(fg_data)]
        
        data_field = str(self.field_combo.currentText())
        if data_field != 'crosstalk':
            data_field = self.analysis[0] + '_' + data_field
        
        if self.analysis[0] == 'pos':
            qc_field = 'ex_qc_pass' if self.analysis[1] == 'ic' else 'in_qc_pass'
        elif self.analysis[0] == 'neg':
            qc_field = 'in_qc_pass' if self.analysis[1] == 'ic' else 'ex_qc_pass'
        fg_x = fg_data[data_field]
        bg_x = bg_data[data_field]
        fg_qc = fg_data[qc_field].copy()
        bg_qc = bg_data[qc_field].copy()

        # record data for later use
        self.fg_x = fg_x
        self.bg_x = bg_x
        self.fg_qc = fg_qc
        self.bg_qc = bg_qc
        self.qc_field = qc_field
        self.fg_data = fg_data
        self.bg_data = bg_data

        # select colors
        pass_brush = pg.mkBrush((255, 255, 255, 80))
        fail_brush = pg.mkBrush((150, 0, 0, 80))
        fg_color = [(pass_brush if qc_pass else fail_brush) for qc_pass in fg_qc]
        bg_color = [(pass_brush if qc_pass else fail_brush) for qc_pass in bg_qc]
        for i in range(len(fg_data)):
            # QC failures are colored red
            if not fg_qc[i]:
                continue

            pulse_n = fg_data['pulse_number'][i] - 1

            # If a pulse number is deselected, then we just mark it as qc-failed and color the point orange
            if not self.pulse_checks[pulse_n].isChecked():
                fg_color[i] = pg.mkBrush(255, 150, 0, 80)
                fg_qc[i] = False
                continue

            # Otherwise, we can color by pulse number if requested
            if self.color_by_pulse_check.isChecked():
                g = pulse_n * 255/7.
                b = 255 - g
                if pulse_n > 7:
                    color = pg.mkColor(0, 0, 0, 0)
                else:
                    color = pg.mkColor(0, g, b, 80)
                fg_color[i] = pg.mkBrush(color)

        # clear old plots
        self.fg_scatter.setData([])
        self.bg_scatter.setData([])
        self.hist_plot.clear()
        self.hist_plot.setTitle("")
        self.clear_trace_plots()
        # clear click-selected plots
        self._clicked_fg_ids = []
        self._clicked_bg_ids = []

        if len(fg_x) == 0:
            return

        # scatter plots of fg/bg data
        fg_y = np.linspace(1.8, 1, len(fg_x))
        bg_y = np.linspace(0.8, 0, len(bg_x))
        self.fg_scatter.setData(fg_x, fg_y, data=fg_data, brush=fg_color)
        self.bg_scatter.setData(bg_x, bg_y, data=bg_data, brush=bg_color)

        # show only qc-passed data in histogram
        fg_x_qc = fg_x[fg_qc] 
        bg_x_qc = bg_x[bg_qc]
        if len(fg_x_qc) == 0 or len(bg_x_qc) == 0:
            return
        
        # plot histograms        
        n_bins = max(5, int(len(fg_x_qc)**0.5))
        all_x = np.concatenate([fg_x_qc, bg_x_qc])
        bins = np.linspace(np.percentile(all_x, 2), np.percentile(all_x, 98), n_bins+1)
        fg_hist = np.histogram(fg_x_qc, bins=bins)
        bg_hist = np.histogram(bg_x_qc, bins=bins)
        self.hist_plot.plot(bg_hist[1], bg_hist[0], stepMode=True, fillLevel=0, brush=(255, 0, 0, 100), pen=0.5)
        self.hist_plot.plot(fg_hist[1], fg_hist[0], stepMode=True, fillLevel=0, brush=(0, 255, 255, 100), pen=1.0)

        # measure distance between distributions
        ks = scipy.stats.ks_2samp(fg_x_qc, bg_x_qc)
        self.hist_plot.setTitle("%s: KS=%0.02g" % (' '.join(self.analysis), ks.pvalue))

        # set event region, but don't update (it's too expensive)
        try:
            self.event_region.sigRegionChangeFinished.disconnect(self.region_changed)
            self.event_region.setRegion([fg_x_qc.min(), fg_x_qc.max()])
        finally:
            self.event_region.sigRegionChangeFinished.connect(self.region_changed)

    def clear_trace_plots(self):
        self.fg_trace_plot.clear()
        self.bg_trace_plot.clear()
        self.selected_fg_traces = []
        self.selected_bg_traces = []

    def region_changed(self):
        """Event selection region changed; replot selected traces
        """
        self.clear_trace_plots()
        rgn = self.event_region.getRegion()

        # Find selected events
        mask = (self.fg_x > rgn[0]) & (self.fg_x < rgn[1])
        fg_ids = self.fg_data['id'][mask]
        self._selected_fg_ids = fg_ids

        mask = (self.bg_x > rgn[0]) & (self.bg_x < rgn[1])
        bg_ids = self.bg_data['id'][mask]
        self._selected_bg_ids = bg_ids

        # generate trace queries
        self.plot_prd_ids(fg_ids, 'fg', trace_list=self.selected_fg_traces, avg=True)
        self.plot_prd_ids(bg_ids, 'bg', trace_list=self.selected_bg_traces, avg=True)

    def replot_all(self):
        self.plot_prd_ids(self._selected_fg_ids, 'fg', pen=None, trace_list=self.selected_fg_traces, avg=True)
        self.plot_prd_ids(self._selected_bg_ids, 'bg', pen=None, trace_list=self.selected_bg_traces, avg=True)
        self.plot_prd_ids(self._clicked_fg_ids, 'fg', pen='y', trace_list=self.clicked_fg_traces)
        self.plot_prd_ids(self._clicked_bg_ids, 'bg', pen='y', trace_list=self.clicked_bg_traces)

    def get_pulse_recs(self, ids, source):
        ids = list(map(int, ids))
        if source == 'fg':
            q = response_query(self.session)
            q = q.join(PulseResponseStrength)
            q = q.filter(PulseResponseStrength.id.in_(ids))
            q = q.add_column(db.PulseResponse.start_time)
            traces = self.selected_fg_traces
            plot = self.fg_trace_plot
        else:
            q = baseline_query(self.session)
            q = q.join(BaselineResponseStrength)
            q = q.filter(BaselineResponseStrength.id.in_(ids))
            q = q.add_column(db.Baseline.start_time)
            traces = self.selected_bg_traces
            plot = self.bg_trace_plot
        
        q = q.join(db.SyncRec).add_column(db.SyncRec.ext_id.label('sync_rec_ext_id'))
        recs = q.all()
        return recs

    def plot_prd_ids(self, ids, source, pen=None, trace_list=None, avg=False):
        """Plot raw or decolvolved PulseResponse data, given IDs of records in
        a PulseResponseStrength table.
        """
        with pg.BusyCursor():
            recs = self.get_pulse_recs(ids, source)
            if len(recs) == 0:
                return

            if source == 'fg':
                traces = self.selected_fg_traces
                plot = self.fg_trace_plot
            else:
                traces = self.selected_bg_traces
                plot = self.bg_trace_plot

            for i in trace_list[:]:
                plot.removeItem(i)
                trace_list.remove(i)
                
            if pen is None:
                alpha = np.clip(1000 / len(recs), 30, 255)
                pen = (255, 255, 255, alpha)
                
            traces = []
            spike_times = []
            spike_values = []
            for rec in recs:
                # Filter by QC unless we selected just a single record
                if len(recs) > 0 and getattr(rec, self.qc_field) is False:
                    continue

                s = {'fg': 'pulse_response', 'bg': 'baseline'}[source]
                result = analyze_response_strength(rec, source=s, lpf=self.lpf_check.isChecked(), 
                                                   remove_artifacts=self.ar_check.isChecked(), bsub=self.bsub_check.isChecked())

                if self.deconv_check.isChecked():
                    trace = result['dec_trace']
                else:
                    trace = result['raw_trace']
                    if self.bsub_check.isChecked():
                        trace = trace - np.median(trace.time_slice(0, 9e-3).data)
                    if self.lpf_check.isChecked():
                        trace = filter.bessel_filter(trace, 500)
                
                spike_values.append(trace.value_at([result['spike_time']])[0])
                if self.align_check.isChecked():
                    trace.t0 = -result['spike_time']
                    spike_times.append(0)
                else:
                    spike_times.append(result['spike_time'])

                traces.append(trace)
                trace_list.append(plot.plot(trace.time_values, trace.data, pen=pen))

            if avg:
                mean = TraceList(traces).mean()
                trace_list.append(plot.plot(mean.time_values, mean.data, pen='g'))
                trace_list[-1].setZValue(10)

            spike_scatter = pg.ScatterPlotItem(spike_times, spike_values, size=4, pen=None, brush=(200, 200, 0))
            spike_scatter.setZValue(-100)
            plot.addItem(spike_scatter)
            trace_list.append(spike_scatter)


def query_all_pairs(classifier=None):
    columns = [
        "connection_strength.*",
        "experiment.id as experiment_id",
        "experiment.acq_timestamp as acq_timestamp",
        "experiment.rig_name",
        "experiment.acsf",
        "slice.species as donor_species",
        "slice.genotype as donor_genotype",
        "slice.age as donor_age",
        "slice.sex as donor_sex",
        "slice.quality as slice_quality",
        "slice.weight as donor_weight",
        "slice.slice_time",
        "pre_cell.ext_id as pre_cell_id",
        "pre_cell.cre_type as pre_cre_type",
        "pre_cell.target_layer as pre_target_layer",
        "pre_morphology.pyramidal as pre_pyramidal",
        "post_cell.ext_id as post_cell_id",
        "post_cell.cre_type as post_cre_type",
        "post_cell.target_layer as post_target_layer",
        "post_morphology.pyramidal as post_pyramidal",
        "pair.synapse",
        "pair.distance",
        "pair.crosstalk_artifact",
        "abs(post_cell.ext_id - pre_cell.ext_id) as electrode_distance",
    ]
    # columns.extend([
    #     "detection_limit.minimum_amplitude",
    # ])

    joins = [
        "join pair on connection_strength.pair_id=pair.id",
        "join cell pre_cell on pair.pre_cell_id=pre_cell.id",
        "join cell post_cell on pair.post_cell_id=post_cell.id",
        "join morphology pre_morphology on pre_morphology.cell_id=pre_cell.id",
        "join morphology post_morphology on post_morphology.cell_id=post_cell.id",
        "join experiment on pair.experiment_id=experiment.id",
        "join slice on experiment.slice_id=slice.id",
    ]
    # joins.extend([
    #     "left join detection_limit on detection_limit.pair_id=pair.id",
    # ])


    query = ("""
    select 
    {columns}
    from connection_strength
    {joins}
    order by acq_timestamp
    """).format(
        columns=", ".join(columns), 
        joins=" ".join(joins),
    )

    session = db.Session()
    df = pandas.read_sql(query, session.bind)

    recs = df.to_records()

    if classifier is None:
        return recs

    # Fit classifier and add results of classifier prediction in to records
    classifier.fit(recs)
    prediction = classifier.predict(recs)
    recs = join_struct_arrays([recs, prediction])
    return recs


pair_classifier = None
def get_pair_classifier(**kwds):
    global pair_classifier
    if pair_classifier is None:
        pair_classifier = PairClassifier(**kwds)
    return pair_classifier
    

def datetime_to_timestamp(d):
    return time.mktime(d.timetuple()) + d.microsecond * 1e-6


def prs_qc():
    """Convenience function for debugging QC: returns (recording, window) arguments
    used for pulse response QC
    """
    global selected_response
    sr = selected_response
    if isinstance(sr, PulseResponseStrength):
        resp = sr.pulse_response
    else:
        resp = sr.baseline
    rec = resp.recording
    start = resp.start_time
    expt = rec.sync_rec.experiment
    nwb = expt.data
    rrec = nwb.contents[rec.sync_rec.ext_id][rec.electrode.device_id]
    dt = rrec['primary'].dt
    dur = len(resp.data) / db.default_sample_rate
    w = [int(start/dt), int((start+dur)/dt)]

    stdev = rrec[w[0]:w[1]]['primary'].std()
    if rrec.clamp_mode == 'ic':
        median = rrec[w[0]:w[1]]['primary'].median()
    else:
        median = rrec[w[0]:w[1]]['command'].median()
    print("Selected:  sweep %d @ %d ms   stdev=%0.3g, med=%0.3g" % (rec.sync_rec.ext_id, np.round(start*1000), stdev, median))

    return rrec, w


class PairClassifier(object):
    """Supervised classifier used to predict whether a cell pair is synaptically connected.

    Input records should be similar to those generated by query_all_pairs()
    """
    def __init__(self, seed=None, use_vc_features=True):
        ic_features = [
            # 'ic_amp_mean',
            # 'ic_amp_stdev',
            # 'ic_amp_ks2samp',  # hurts perfornamce
            'ic_deconv_amp_mean',
            # 'ic_deconv_amp_stdev',
            'ic_deconv_amp_ks2samp',
            'ic_latency_mean',
            'ic_latency_stdev',
            'ic_latency_ks2samp',
            # 'ic_crosstalk_mean',
            'ic_fit_amp',
            'ic_fit_xoffset',
            'ic_fit_yoffset',
            'ic_fit_rise_time',
            #'ic_fit_rise_power',
            'ic_fit_decay_tau',
            #'ic_fit_exp_amp',
            'ic_fit_nrmse',

        ]
        vc_features = [
            # 'vc_amp_mean',  # hurts perfornamce
            # 'vc_amp_stdev',
            'vc_amp_ks2samp',
            # 'vc_deconv_amp_mean',
            # 'vc_deconv_amp_stdev',
            # 'vc_deconv_amp_ks2samp',  # hurts performance
            # 'vc_latency_mean',
            # 'vc_latency_stdev',  # hurts performance
            'vc_latency_ks2samp',

            'vc_fit_amp',
            'vc_fit_xoffset',
            'vc_fit_yoffset',
            'vc_fit_rise_time',
            #'vc_fit_rise_power',
            'vc_fit_decay_tau',
            #'vc_fit_exp_amp',
            'vc_fit_nrmse',
        ]
        general_features = [
            # 'electrode_distance',
        ]

        self.features = general_features + ic_features
        if use_vc_features:
             self.features.extend(vc_features)

        # Random seed used when shuffling training/test inputs
        self.seed = seed

    def fit(self, recs):
        ids = recs['id']

        # Select features from records
        features = np.array([tuple(r) for r in recs[self.features]])

        # QC bad records; don't train on these
        mask = (recs['ic_n_samples'] > 100) & (recs['ic_crosstalk_mean'] < 60e-6)

        x = features[mask]
        y = recs['synapse'][mask].astype(bool)
        ids = ids[mask]

        # shuffle
        order = np.arange(len(y))
        rand = np.random.RandomState(seed=self.seed)
        rand.shuffle(order)
        x = x[order]
        y = y[order]
        ids = ids[order]

        # mask out nan/inf records
        mask2 = np.all(np.isfinite(x), axis=1) & np.isfinite(y)
        x = x[mask2]
        y = y[mask2]
        ids = ids[mask2]
        mask[mask] = mask2

        # prescale records, keep the scaler for later processing
        scaler = sklearn.preprocessing.StandardScaler().fit(x)
        x = scaler.transform(x)
        self.scaler = scaler

        # split into training and test sets
        # select training set from connected and non-connected separately to ensure 
        # we get enough connected examples in the training set
        train_mask = np.zeros(len(y), dtype='bool')
        syn = np.argwhere(y)
        n = len(syn) // 4 * 3
        train_mask[syn[:n]] = True
        other = np.arange(len(y))[~train_mask]
        rand.shuffle(other)
        # 1/5 of training data is connected, the rest is unconnected
        train_mask[other[:n*5]] = True

        train_x = x[train_mask]
        train_y = y[train_mask]
        test_x = x[~train_mask]
        test_y = y[~train_mask]
        print("Train: %d  test: %d   random seed: %s" % (len(train_y), len(test_y), self.seed))

        # build and fit the classifier
        #clf = sklearn.svm.LinearSVC()
        clf = sklearn.svm.SVC(C=1, class_weight='balanced', coef0=0.0,
        decision_function_shape='ovr', degree=3, gamma='auto', kernel='rbf',
        max_iter=-1, probability=True, random_state=self.seed)
        # clf = sklearn.ensemble.RandomForestClassifier()

        hyper_params = [{'C': [1, 10, 100, 1000], 'gamma': [0.1, 0.01, 0.001, 0.0001]}]
        clf = sklearn.model_selection.GridSearchCV(clf, hyper_params)
        self.clf = clf

        clf.fit(train_x, train_y)
        print("Classifier best hyperparameters:", clf.best_params_)

        def test(clf, test_x, test_y):
            pred = clf.predict(test_x)
            def pr(name, pred, test_y):
                print("%s:  %d/%d  %0.2f%%" % (name, (pred & test_y).sum(), test_y.sum(), 100 * (pred & test_y).sum() / test_y.sum()))
                
            pr(" True positive", pred, test_y)
            pr("False positive", pred, ~test_y)
            pr(" True negative", ~pred, ~test_y)
            pr("False negative", ~pred, test_y)

        print("TRAIN:")
        test(clf, train_x, train_y)

        print("TEST:")
        test(clf, test_x, test_y)

        # calculate approximate probability threshold between classes
        pred = self.predict(recs)
        min_prob = pred[pred['prediction'] == True]['confidence'].min()
        max_prob = pred[pred['prediction'] == False]['confidence'].max()
        self.prob_threshold = (min_prob + max_prob) / 2.

    def predict(self, recs=None):
        """Predict connectivity for a sequence of records output from analyze_response_strength

        Input may be a structured array or list of dicts.
        """
        # Select features from records
        if isinstance(recs, np.ndarray):
            features = np.array([tuple(r) for r in recs[self.features]])
        else:
            features = np.array([tuple(map(r.__getitem__, self.features)) for r in recs])

        # prepare ouptut array
        result = np.empty(len(features), dtype=[('prediction', float), ('confidence', float)])
        result[:] = np.nan

        # mask out inf/nan records
        mask = np.all(np.isfinite(features), axis=1)

        # scale masked records
        norm_features = self.scaler.transform(features[mask])

        # fill in output array
        result['prediction'][mask] = self.clf.predict(norm_features)
        # result['confidence'][mask] = clf.decision_function(norm_features)
        result['confidence'][mask] = self.clf.predict_proba(norm_features)[:,1]
        # result['confidence'][mask] = clf.predict_proba(norm_features)[:,1]
        assert np.isfinite(result[mask]['confidence']).sum() > 0
        return result


def join_struct_arrays(arrays):
    """Join two structured arrays together.

    This is the most inefficient possible approach, but other methods
    don't work well with object dtypes.
    """
    dtype = []
    for arr in arrays:
        for name in arr.dtype.names:
            dtype.append((str(name), arr.dtype.fields[name][0].str))
    arr = np.empty(len(arrays[0]), dtype=dtype)
    for i in range(len(arr)):
        v = ()
        for a in arrays:
            v = v + tuple(a[i])
        arr[i] = v
    return arr


class PairScatterPlot(pg.QtCore.QObject):
    """Create a ScatterPlotWidget that displays results selected with query_all_pairs()
    """

    pair_clicked = pg.QtCore.Signal(object)

    def __init__(self, recs):
        pg.QtCore.QObject.__init__(self)

        # Load all records into scatter plot widget
        spw = pg.ScatterPlotWidget()
        spw.style['symbolPen'] = None
        
        spw.resize(1000, 800)
        spw.show()

        fields = [
            ('synapse', {'mode': 'enum', 'values': [True, False, None]}),
            ('synapse_type', {'mode': 'enum', 'values': ['in', 'ex']}),
            ('prediction', {'mode': 'enum', 'values': [True, False, None]}),
            ('confidence', {}),
            ('pre_cre_type', {'mode': 'enum', 'values': list(set(recs['pre_cre_type']))}),
            ('post_cre_type', {'mode': 'enum', 'values': list(set(recs['post_cre_type']))}),
            ('pre_pyramidal', {'mode': 'enum', 'values': list(set(recs['pre_pyramidal']))}),
            ('post_pyramidal', {'mode': 'enum', 'values': list(set(recs['post_pyramidal']))}),
            ('pre_target_layer', {'mode': 'enum'}),
            ('post_target_layer', {'mode': 'enum'}),
            ('ic_n_samples', {}),
            ('vc_n_samples', {}),
            ('rig_name', {'mode': 'enum', 'values': list(set(recs['rig_name']))}),
            ('acsf', {'mode': 'enum', 'values': list(set(recs['acsf']))}),
            ('acq_timestamp', {}),
            ('crosstalk_artifact', {'units': 'V'}),
            ('electrode_distance', {}),
            ('slice_quality', {'mode': 'enum', 'values': list(range(1,6))}),
        ]
        fnames = [f[0] for f in fields]
        for f in recs.dtype.names:
            if f in fnames:
                continue
            if 'amp' in f:
                if f[:2] == 'vc':
                    fields.append((f, {'units': 'A'}))
                else:
                    fields.append((f, {'units': 'V'}))
            elif 'latency' in f:
                fields.append((f, {'units': 's'}))
            else:
                fields.append((f, {}))
                
        spw.setFields(fields)
        spw.setData(recs)

        spw.sigScatterPlotClicked.connect(self._conn_clicked)

        # Set up scatter plot widget defaults
        spw.setSelectedFields('ic_base_deconv_amp_mean', 'ic_deconv_amp_mean')
        spw.filter.addNew('synapse')
        ch = spw.filter.addNew('ic_crosstalk_mean')
        ch['Min'] = -1
        ch['Max'] = 60e-6
        ch = spw.filter.addNew('rig_name')
        ch = spw.colorMap.addNew('synapse')
        ch['Values', 'True'] = pg.mkColor('y')
        ch['Values', 'False'] = pg.mkColor(200, 200, 200, 150)
        ch['Values', 'None'] = pg.mkColor(200, 100, 100)
        ch = spw.colorMap.addNew('ic_latency_mean')
        ch['Operation'] = 'Add'
        ch['Max'] = 3e-3
        cm = pg.ColorMap([0, 1], [[0, 0, 255, 255], [0, 0, 0, 255]])
        ch.setValue(cm)

        self.spw = spw

    def _conn_clicked(self, spw, points):
        self.spw.setSelectedPoints([points[0]])
        d = points[0].data()
        self.selected = d
        self.pair_clicked.emit(d['pair_id'])


class PairView(pg.QtCore.QObject):
    """For browsing and analyzing pairs. 

    Contains an ExperimentBrowser and a number of analyzers showing response properties of a selected pair
    """
    def __init__(self):
        pg.QtCore.QObject.__init__(self)

        # global session for querying from DB
        self.session = db.Session()

        win = pg.QtGui.QSplitter()
        win.setOrientation(pg.QtCore.Qt.Horizontal)
        win.resize(1000, 800)
        win.show()
        
        b = ExperimentBrowser()
        win.addWidget(b)
        
        rs_plots = ResponseStrengthPlots(self.session)
        win.addWidget(rs_plots)

        b.itemSelectionChanged.connect(self._selected)            
        b.doubleClicked.connect(self._dbl_clicked)

        self.win = win
        self.rs_plots = rs_plots
        self.browser = b
        self.nwb_viewer = MultipatchNwbViewer()

    def select_pair(self, pair_id):
        self.browser.select_pair(pair_id)

    def _dbl_clicked(self, index):
        """Item double clicked; load in NWB viewer
        """
        with pg.BusyCursor():
            item = self.browser.itemFromIndex(index)[0]
            self.nwb_viewer.load_nwb(item.expt.nwb_cache_file)
            self.nwb_viewer.show()

    def _selected(self, *args):
        """A pair was selected; update the event plots
        """
        sel = self.browser.selectedItems()
        self.selected = sel
        if len(sel) == 0:
            return
        sel = sel[0]
        if hasattr(sel, 'pair'):
            pair = sel.pair
            expt = pair.experiment
            self.rs_plots.load_conn(pair)

        sec = sel.expt.acq_timestamp
        src = sel.expt.source_experiment
        print("======================================")
        if src.entry is not None:
            src.entry.print_tree()
        else:
            print(open(src.source_id[0], 'rb').read())
        print("------------------------------")
        print("Original path:", src)
        print("Server path:", sel.expt.original_path)
        if hasattr(sel, 'pair'):
            print("ID: %.3f  %d->%d" % (sec, pair.pre_cell.ext_id, pair.post_cell.ext_id))
            conn = pair.connection_strength
            cls = get_pair_classifier()
            f = {k: getattr(conn, k) for k in cls.features}
            print(f)
            print(cls.predict([f]))
        else:
            print("ID: %s" % sec)
        



def str_analysis_result_table(results, recs):
    """Convert output of strength_analysis.analyze_response_strength to look like
    the result was queried from the DB using get_amps() or get_baseline()
    """
    dtype = [
        ('id', int),
        ('pos_amp', float),
        ('neg_amp', float),
        ('pos_dec_amp', float),
        ('neg_dec_amp', float),
        ('pos_dec_latency', float),
        ('neg_dec_latency', float),
        ('crosstalk', float),
        ('ex_qc_pass', bool),
        ('in_qc_pass', bool),
        ('clamp_mode', object),
        ('pulse_number', int),
        ('max_dvdt_time', float),
        ('response_start_time', float),
        ('data', object),
        ('rec_start_time', float),
    ]
    
    table = np.empty(len(recs), dtype=dtype)
    for i,rec in enumerate(recs):
        for key in ['ex_qc_pass', 'in_qc_pass', 'clamp_mode', 'data']:
            table[i][key] = getattr(rec, key)
        result = results[i]
        for key,val in result.items():
            if key in table.dtype.names:
                table[i][key] = val
        table[i]['max_dvdt_time'] = 10e-3
        table[i]['response_start_time'] = 0
    return table


class RecordWrapper(object):
    """Wraps records returned from DB so that we can override some values.
    """
    def __init__(self, rec):
        self._rec = rec
    def __getattr__(self, attr):
        return getattr(self._rec, attr)
        

def simulate_response(fg_recs, bg_results, amp, rtime, seed=None):
    if seed is not None:
        np.random.seed(seed)

    dt = 1.0 / db.default_sample_rate
    t = np.arange(0, 15e-3, dt)
    template = Psp.psp_func(t, xoffset=0, yoffset=0, rise_time=rtime, decay_tau=15e-3, amp=1, rise_power=2)

    r_amps = scipy.stats.binom.rvs(p=0.2, n=24, size=len(fg_recs)) * scipy.stats.norm.rvs(scale=0.3, loc=1, size=len(fg_recs))
    r_amps *= amp / r_amps.mean()
    r_latency = np.random.normal(size=len(fg_recs), scale=200e-6, loc=13e-3)
    fg_results = []
    traces = []
    fg_recs = [RecordWrapper(rec) for rec in fg_recs]  # can't modify fg_recs, so we wrap records with a mutable shell
    for k,rec in enumerate(fg_recs):
        rec.data = rec.data.copy()
        start = int(r_latency[k] * db.default_sample_rate)
        length = len(rec.data) - start
        rec.data[start:] += template[:length] * r_amps[k]

        fg_result = analyze_response_strength(rec, 'baseline')
        fg_results.append(fg_result)

        traces.append(Trace(rec.data, sample_rate=db.default_sample_rate))
        traces[-1].amp = r_amps[k]
    fg_results = str_analysis_result_table(fg_results, fg_recs)
    conn_result = analyze_pair_connectivity({('ic', 'fg'): fg_results, ('ic', 'bg'): bg_results, ('vc', 'fg'): [], ('vc', 'bg'): []}, sign=1)
    return conn_result, traces


def simulate_connection(fg_recs, bg_results, classifier, amp, rtime, n_trials=8):
    """Run repeated simulation trials adding a synthetic PSP to recorded background noise.
    """
    import pyqtgraph.multiprocess as mp
    result = {'results': [], 'rise_time': rtime, 'amp': amp}

    sim_results = [None] * n_trials
    with mp.Parallelize(range(n_trials), results=sim_results, workers=8) as tasker:
        for ii in tasker:
            tasker.results[ii] = simulate_response(fg_recs, bg_results, amp, rtime, seed=ii)

    for k in range(len(sim_results)):
        conn_result, traces = sim_results[k]

        result['results'].append(conn_result)
        result['traces'] = traces
        # print(conn_result)
        # print(dict([(k, conn_result[k]) for k in classifier.features]))

    pred = classifier.predict(result['results'])
    result['predictions'] = pred['prediction']
    result['confidence'] = pred['confidence']
    # print("\nrise time:", rtime, " amplitude:", amp)
    # print(pred)
        
    return result



if __name__ == '__main__':
    import user

    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=0, help="Seed used to randomize classifier inputs")
    args = parser.parse_args(sys.argv[1:])

    pg.dbg()

    # Load records on all pairs and train a classifier to predict connections
    classifier = get_pair_classifier(seed=None if args.seed < 0 else args.seed)
    recs = query_all_pairs(classifier)

    # show all records in scatter plot
    spw = PairScatterPlot(recs)

    # add another window for analyzing selected pairs
    pair_view = PairView()
    spw.pair_clicked.connect(pair_view.select_pair)


    # Print a little report about synapses that may have been misclassified
    session = db.Session()

    fn_mask = (recs['confidence'] > 0.2) & (recs['synapse'] == False)
    fp_mask = (recs['confidence'] < 0.2) & (recs['synapse'] == True)
    for mask, name in [(fn_mask, 'negatives'), (fp_mask, 'positives')]:
        print("\n================ Possible false %s: ==============\n" % name)
        for rec in recs[mask]:
            pid = rec['pair_id']
            pair = session.query(db.Pair).filter(db.Pair.id==pid).all()[0]
            pre_cell = pair.pre_cell
            post_cell = pair.post_cell

            # just select excitatory for now
            if pre_cell.cre_type not in constants.EXCITATORY_CRE_TYPES and pre_cell.morphology.pyramidal is not True:
                continue
            if post_cell.cre_type not in constants.EXCITATORY_CRE_TYPES and post_cell.morphology.pyramidal is not True:
                continue
            
            print("{:s} {:0.3f} {:d} {:d} {:15s} {:20s}  (L{:<3s} {:7s} {:6s}) (L{:<3s} {:7s} {:6s})".format(
                pair.experiment.rig_name, pair.experiment.acq_timestamp, pre_cell.ext_id, post_cell.ext_id, pair.experiment.internal, pair.experiment.acsf,
                pre_cell.target_layer, pre_cell.cre_type, {True: 'pyr', None: '?', False: 'nonpyr'}[pre_cell.morphology.pyramidal], 
                post_cell.target_layer, post_cell.cre_type, {True: 'pyr', None: '?', False: 'nonpyr'}[post_cell.morphology.pyramidal]
            ))


    if sys.flags.interactive == 0:
        pg.mkQApp().exec_()