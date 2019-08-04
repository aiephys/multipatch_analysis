import pyqtgraph as pg


class PairAvgFitUi(object):
    def __init__(self):
        self.widget = pg.QtGui.QSplitter(pg.QtCore.Qt.Horizontal)
        self.plots = {}
        for clamp_mode in ('ic', 'vc'):
            units = {'ic': 'V', 'vc': 'A'}[clamp_mode]
            for holding in (-70, -55):
                self.plots[clamp_mode, holding] = PulseResponsePlot("%s %d" % (clamp_mode, holding), units)

    def show_pulse_responses(self, prs):
        for (clamp_mode, holding), plot in self.plots.items():
            plot.show_pulse_responses(prs[clamp_mode, holding]['qc_pass'], prs[clamp_mode, holding]['qc_fail'])

    def show_data_notes(self, notes):
        for (clamp_mode, holding), plot in self.plots.items():
            plot.show_expected_fit(notes[clamp_mode, holding])

    def show_fit_results(self, results, average):



class PulseResponsePlot(object):
    def __init__(self, title, units):
        self.grid = PlotGrid()
        self.grid.set_shape(2, 1)
        self.grid.grid.ci.layout.setRowStretchFactor(0, 3)
        self.grid.grid.ci.layout.setRowStretchFactor(1, 8)
        self.spike_plot = self.grid[0, 0]
        self.response_plot = self.grid[1, 0]

        self.spike_plot.hideAxis('bottom')
        self.spike_plot.setTitle(title)
        self.spike_plot.setLabel('left', text="presynaptic spike")
        self.spike_plot.addLine(x=0)

        self.response_plot.setLabel('bottom', text='Time from spike', units='s')
        self.response_plot.setXLink(self.spike_plot)
        # self.response_plot.setLabel('left', text="%d holding" % int(holding), units=units)
        
        self.items = []

    def plot_responses(self, pulse_responses):
        self.plot_traces(pulse_responses)
        self.plot_spikes(pulse_responses)

    def plot_traces(self, pulse_responses):  
        for i, holding in enumerate(pulse_responses.keys()):
            for qc, prs in pulse_responses[holding].items():
                if len(prs) == 0:
                    continue
                traces = []
                for pr in prs:
                    trace = pr.post_tseries
                    item = self.trace_plots[i].plot(trace.time_values, trace.data, pen=self.qc_color[qc])
                    if qc == 'qc_fail':
                        item.setZValue(-10)
                    self.items.append(item)
                    traces.append(trace)
                if qc == 'qc_pass':
                    grand_trace = TSeriesList(traces).mean()
                    item = self.trace_plots[i].plot(grand_trace.time_values, grand_trace.data, pen={'color': 'b', 'width': 2})
                    self.items.append(item)
            self.trace_plots[i].autoRange()
            self.trace_plots[i].setXRange(-5e-3, 10e-3)
            # y_range = [grand_trace.data.min(), grand_trace.data.max()]
            # self.plots[i].setYRange(y_range[0], y_range[1], padding=1)

    def plot_spikes(self, pulse_responses):
        for i, holding in enumerate(pulse_responses.keys()):
            for qc, prs in pulse_responses[holding].items():
                if len(prs) == 0:
                    continue
                for pr in prs:
                    spike = pr.pre_tseries
                    item = self.spike_plots[i].plot(spike.time_values, spike.data, pen=self.qc_color[qc])
                    if qc == 'qc_fail':
                        item.setZValue(-10)
                    self.items.append(item)

    def plot_fit(self, trace, fit, holding, fit_pass=False):
        if holding == '-55':
            if self.fit_item_55 is not None:
                self.trace_plots[0].removeItem(self.fit_item_55)
            self.fit_item_55 = pg.PlotDataItem(trace.time_values, fit, name='-55 holding', pen={'color': self.fit_color[fit_pass], 'width': 3})
            self.trace_plots[0].addItem(self.fit_item_55)
        
        elif holding == '-70':
            if self.fit_item_70 is not None:
                self.trace_plots[1].removeItem(self.fit_item_70)
            self.fit_item_70 = pg.PlotDataItem(trace.time_values, fit, name='-70 holding', pen={'color': self.fit_color[fit_pass], 'width': 3})
            self.trace_plots[1].addItem(self.fit_item_70)

    def color_fit(self, name, value):
        if '-55' in name:
            if self.fit_item_55 is not None:
                self.fit_item_55.setPen({'color': self.fit_color[value], 'width': 3})
        if '-70' in name:
            if self.fit_item_70 is not None:
                self.fit_item_70.setPen({'color': self.fit_color[value], 'width': 3})

    def clear_plots(self):
        for item in self.items + [self.fit_item_55, self.fit_item_70]:
            if item is None:
                continue
            item.scene().removeItem(item)
        self.items = []
        
        self.plots[-1].autoRange()
        self.plots[-1].setXRange(-5e-3, 10e-3)
        self.fit_item_70 = None
        self.fit_item_55 = None
        
