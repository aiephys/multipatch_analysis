from __future__ import division, print_function

import time
import pyqtgraph as pg
from pyqtgraph import parametertree
from pyqtgraph.parametertree import Parameter
from neuroanalysis.ui.user_test import UserTestUi

class PulseResponseQCUI(object):
    def __init__(self, title=None):
        self.pw = pg.GraphicsLayoutWidget()
        self.in_plt = self.pw.addPlot(title='Inhibitory Pulse QC')
        self.ex_plt = self.pw.addPlot(row=1, col=0, title='Excitatory Pulse QC')
        self.qc_text = pg.DataTreeWidget()
        self.qc_text.setHeaderLabels(['QC failures', 'type', 'value'])
        self.console = pg.console.ConsoleWidget()

        self.widget = pg.QtGui.QSplitter(pg.QtCore.Qt.Vertical)
        self.widget.addWidget(self.pw)        
        self.widget.addWidget(self.qc_text)
        self.widget.addWidget(self.console)
        self.widget.resize(1000, 900)
        self.widget.setSizes([800, 200])
        self.widget.show()

    def clear(self):
        self.in_plt.clear()
        self.ex_plt.clear()

    def show_result(self, pulse_response, clamp_mode, ex_pass, in_pass, failures):
        pr = pulse_response['primary']
        base = pr.median()
        base_std = pr.std()
        pen = pg.mkPen('y', style=pg.QtCore.Qt.DotLine)
        for plt, qc_pass in zip([self.in_plt, self.ex_plt], [in_pass, ex_pass]):
            trace_color = 'g' if qc_pass else 'r'
            plt.plot(pr.time_values, pr.data, pen=trace_color)
            if clamp_mode == 'vc':
                plt.addLine(y=base + base_std, pen=pen)
                plt.addLine(y=base - base_std, pen=pen)
                plt.addLine(y=base + 15e-12, pen='y')
                plt.addLine(y=base - 15e-12, pen='y')
                plt.setLabel('left', text='Voltage Clamp', units='A')
            if clamp_mode == 'ic':
                plt.addLine(y=base + base_std, pen=pen)
                plt.addLine(y=base - base_std, pen=pen)
                plt.addLine(y=base + 1.5e-3, pen='y')
                plt.addLine(y=base - 1.5e-3, pen='y')
                plt.addLine(y=-40e-3, pen='m')
                plt.addLine(y=-45e-3, pen='c')
                plt.setLabel('left', text='Current Clamp', units='V')
                plt.setYRange(-80e-3, -60e-3)
        if clamp_mode == 'ic':
            self.ex_plt.addLine(y=-85e-3, pen='c')
            self.in_plt.addLine(y=-60e-3, pen='c')

        # qc_fail_text = '\n'.join(failures)
        # self.qc_text.setValue(qc_fail_text)
        self.qc_text.setData(failures)


class PulseResponseQCTestUi(UserTestUi):
    """UI for manually pass/failing pulse response QC unit tests.
    """
    def __init__(self):
        expected_display = PulseResponseQCUI('expected result')
        current_display = PulseResponseQCUI('current result')
        UserTestUi.__init__(self, expected_display, current_display)