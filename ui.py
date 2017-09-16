from collections import OrderedDict
import os, sys, subprocess
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.miesnwb import MiesNwb


class SynapseTreeWidget(QtGui.QTreeWidget):
    """Displays a list of known synapses, given a lit of experiments.
    
    Provides ui for filtering and selecting.
    """
    def __init__(self, expts, parent=None):
        QtGui.QTreeWidget.__init__(self, parent)
        self.setColumnCount(3)
        self.expts = expts
        for syn in expts.connection_summary():
            item = SynapseTreeItem(syn['expt'], syn['cells'])
            self.addTopLevelItem(item)

    
class SynapseTreeItem(QtGui.QTreeWidgetItem):
    def __init__(self, expt, cells):
        self.expt = expt
        self.cells = cells

        fields = [
            str(expt.summary_id),
            "%d - %d" % (cells[0].cell_id, cells[1].cell_id),
            "%s - %s" % (cells[0].cre_type, cells[1].cre_type),
            
        ]
        
        QtGui.QTreeWidgetItem.__init__(self, fields)


class ExperimentInfoWidget(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        
        self.expt = None
        
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        
        self.info_tree = pg.DataTreeWidget()
        self.layout.addWidget(self.info_tree, self.layout.rowCount(), 0)
        
        self.biocytin_btn = QtGui.QPushButton('biocytin image...')
        self.biocytin_btn.clicked.connect(self.show_biocytin)
        self.layout.addWidget(self.biocytin_btn, self.layout.rowCount(), 0)
        
    def set_experiment(self, expt):
        self.expt = expt
        
        info = OrderedDict([
            ('date', str(expt.date)),
            ('specimen', expt.specimen_id),
            ('age', expt.age),
        ])
          
        self.info_tree.setData(info)
        
    def show_biocytin(self):
        if sys.platform == 'win32':
            subprocess.Popen([r'C:\Program Files (x86)\Mozilla Firefox\firefox.exe', self.expt.biocytin_image_url])
        else:
            os.system("firefox " + self.expt.biocytin_image_url)


class ExperimentTimeline(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        self.plots = PlotGrid()
        self.layout.addWidget(self.plots, 0, 0)

        self.ptree = pg.parametertree.ParameterTree(showHeader=False)
        self.layout.addWidget(self.ptree, 0, 1)
        self.ptree.setMaximumWidth(250)
        
        self.params = pg.parametertree.Parameter.create(name='params',
            type='group', addText='Add electrode')
        self.params.addNew = self.add_electrode_clicked  # monkey!
        self.ptree.setParameters(self.params, showTop=False)
        
    def add_electrode_clicked(self):
        elec = ElectrodeParameter(self)
        self.params.addChild(elec)
        
    def load_experiment(self, nwb_handle):
        self.nwb_handle = nwb_handle
        self.nwb = MiesNwb(nwb_handle.name())
        
        # load all recordings
        recs = {}
        for srec in self.nwb.contents:
            for chan in srec.devices:
                recs.setdefault(chan, []).append(srec[chan])

        chans = sorted(recs.keys())
        self.plots.set_shape(len(chans), 1)
        self.plots.setXLink(self.plots[0, 0])
        
        # find time of first recording
        start_time = min([rec[0].start_time for rec in recs.values()])
        self.start_time = start_time
        end_time = max([rec[-1].start_time for rec in recs.values()])
        self.plots.setXRange(0, (end_time-start_time).seconds)
        
        # plot all recordings
        for i,chan in enumerate(chans):
            n_recs = len(recs[chan])
            times = np.empty(n_recs)
            i_hold = np.empty(n_recs)
            v_hold = np.empty(n_recs)
            v_noise = np.empty(n_recs)
            i_noise = np.empty(n_recs)
            
            # load QC metrics for all recordings
            for j,rec in enumerate(recs[chan]):
                dt = (rec.start_time - start_time).seconds
                times[j] = dt
                v_hold[j] = rec.baseline_potential
                i_hold[j] = rec.baseline_current
                if rec.clamp_mode == 'vc':
                    v_noise[j] = np.nan
                    i_noise[j] = rec.baseline_rms_noise
                else:
                    v_noise[j] = rec.baseline_rms_noise
                    i_noise[j] = np.nan
                    
            # scale all qc metrics to the range 0-1
            pass_brush = pg.mkBrush('b')
            fail_brush = pg.mkBrush('r')
            v_hold = (v_hold + 60e-3) / 20e-3
            i_hold = i_hold / 400e-12
            v_noise = v_noise / 5e-3
            i_noise = i_noise / 100e-12
            
            plt = self.plots[i, 0]
            plt.setLabels(left=("Ch %02d" % chan))
            for data,symbol in [(np.zeros_like(times), 'o'), (v_hold, 't'), (i_hold, 'x'), (v_noise, 't1'), (i_noise, 'x')]:
                brushes = np.where(np.abs(data) > 1.0, fail_brush, pass_brush)
                plt.plot(times, data, pen=None, symbol=symbol, symbolPen=None, symbolBrush=brushes)

        # automatically select electrode regions
        site_info = self.nwb_handle.parent().info()
        for i in range(1, 9):
            hs_state = site_info.get('Headstage %d'%i, None)
            if hs_state is None:
                continue
            
            start = recs[i-1][0].start_time
            stop = recs[i-1][-1].start_time
            elec = ElectrodeParameter(self, i, start, stop)
            self.params.addChild(elec)
            
        

class ElectrodeParameter(pg.parametertree.parameterTypes.GroupParameter):
    def __init__(self, ui):
        self.ui = ui
        params = [
            {'name': 'channel', 'type': 'list', 'values': []},
            {'name': 'status', 'type': 'list', 'values': ['No seal', 'Low seal', 'GOhm seal', 'Technical failure']},
        ]
        pg.parametertree.parameterTypes.GroupParameter.__init__(self, name="Electrode?", children=params, removable=True)
        
        