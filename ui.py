from collections import OrderedDict
import os, sys, subprocess, datetime
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.miesnwb import MiesNwb
from config import n_headstages


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
            str(expt.uid),
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
        self.channels = None
        
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        self.plots = PlotGrid()
        self.layout.addWidget(self.plots, 0, 0)

        self.ptree = pg.parametertree.ParameterTree(showHeader=False)
        self.layout.addWidget(self.ptree, 0, 1)
        self.ptree.setMaximumWidth(250)
        
        self.params = pg.parametertree.Parameter.create(name='params',
            type='group', addText='Add pipette')
        self.params.addNew = self.add_pipette_clicked  # monkey!
        self.ptree.setParameters(self.params, showTop=False)
        
    def list_channels(self):
        return self.channels
        
    def get_channel_plot(self, chan):
        i = self.channels.index(chan)
        return self.plots[i, 0]
        
    def add_pipette_clicked(self):
        self.add_pipette(channel=self.channels[0], start=0, stop=500)
        
    def remove_pipettes(self):
        for ch in self.params.children():
            self.params.removeChild(ch)
            ch.region.scene().removeItem(ch.region)
        
    def load_experiment(self, nwb_handle):
        self.nwb_handle = nwb_handle
        self.nwb = MiesNwb(nwb_handle.name())
        
        # load all recordings
        recs = {}
        for srec in self.nwb.contents:
            for chan in srec.devices:
                recs.setdefault(chan, []).append(srec[chan])

        chans = sorted(recs.keys())
        self.channels = chans
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
            pass_brush = pg.mkBrush(100, 100, 255, 200)
            fail_brush = pg.mkBrush(255, 0, 0, 200)
            v_hold = (v_hold + 60e-3) / 20e-3
            i_hold = i_hold / 400e-12
            v_noise = v_noise / 5e-3
            i_noise = i_noise / 100e-12
            
            plt = self.plots[i, 0]
            plt.setLabels(left=("Ch %d" % chan))
            for data,symbol in [(np.zeros_like(times), 'o'), (v_hold, 't'), (i_hold, 'x'), (v_noise, 't1'), (i_noise, 'x')]:
                brushes = np.where(np.abs(data) > 1.0, fail_brush, pass_brush)
                plt.plot(times, data, pen=None, symbol=symbol, symbolPen=None, symbolBrush=brushes)

        # automatically select electrode regions
        self.remove_pipettes()
        site_info = self.nwb_handle.parent().info()
        for i in self.channels:
            hs_state = site_info.get('Headstage %d'%i, None)
            if hs_state is None:
                continue
            status = {
                'NS': 'No seal',
                'LS': 'Low seal',
                'GS': 'GOhm seal',
                'TF': 'Technical failure',
            }[hs_state]
            start = (recs[i][0].start_time - start_time).seconds - 1
            stop = (recs[i][-1].start_time - start_time).seconds + 1
            
            # assume if we got more than two recordings, then a cell was present.
            got_cell = len(recs[i]) > 2
            
            self.add_pipette(i, start, stop, status=status, got_cell=got_cell)
            
    def add_pipette(self, channel, start, stop, status=None, got_cell=None):
        elec = PipetteParameter(self, channel, start, stop, status=status, got_cell=got_cell)
        self.params.addChild(elec, autoIncrementName=True)
        elec.child('channel').sigValueChanged.connect(self._pipette_channel_changed)
        elec.region.sigRegionChangeFinished.connect(self._pipette_region_changed)
        self._pipette_channel_changed(elec.child('channel'))
        
    def _pipette_channel_changed(self, param):
        plt = self.get_channel_plot(param.value())
        plt.addItem(param.parent().region)
        self._rename_pipettes()
        
    def _pipette_region_changed(self):
        self._rename_pipettes()
        
    def _rename_pipettes(self):
        # sort electrodes by channel
        elecs = {}
        for elec in self.params.children():
            elecs.setdefault(elec['channel'], []).append(elec)
        
        for chan in elecs:
            # sort all electrodes on this channel by start time
            chan_elecs = sorted(elecs[chan], key=lambda e: e.region.getRegion()[0])
            
            # assign names
            for i,elec in enumerate(chan_elecs):
                # rename all first to avoid name colisions
                elec.setName('rename%d' % i)
            for i,elec in enumerate(chan_elecs):
                # If there are multiple electrodes on this channel, then 
                # each extra electrode increments its name by the number of 
                # headstages (for example, on AD channel 3, the first electrode
                # is called "Electrode 4", and on an 8-headstage system, the
                # second electrode will be "Electrode 12").
                e_id = (chan+1) + (i*n_headstages)
                elec.id = e_id
                elec.setName('Pipette %d' % e_id)

    def save(self):
        state = []
        for elec in self.params.children():
            rgn = elec.region.getRegion()
            start = self.start_time + datetime.timedelta(seconds=rgn[0])
            stop = self.start_time + datetime.timedelta(seconds=rgn[1])
            state.append({
                'id': elec.id,
                'status': elec['status'],
                'got_cell': elec['got cell'],
                'channel': elec['channel'],
                'start': start,
                'stop': stop,
            })
        return state


class PipetteParameter(pg.parametertree.parameterTypes.GroupParameter):
    def __init__(self, ui, channel, start, stop, status=None, got_cell=None):
        self.ui = ui
        params = [
            {'name': 'channel', 'type': 'list', 'values': ui.list_channels()},
            {'name': 'status', 'type': 'list', 'values': ['No seal', 'Low seal', 'GOhm seal', 'Technical failure']},
            {'name': 'got cell', 'type': 'bool'},
        ]
        pg.parametertree.parameterTypes.GroupParameter.__init__(self, name="Pipette?", children=params, removable=True)
        self.child('got cell').sigValueChanged.connect(self._got_cell_changed)
        
        region = [0, 500]
        if start is not None:
            region[0] = start
        if stop is not None:
            region[1] = stop

        self.region = pg.LinearRegionItem(region)
        self.region.setZValue(-10)

        if channel is not None:
            self['channel'] = channel
        if status is not None:
            self['status'] = status
        if got_cell is not None:
            self['got cell'] = got_cell
            
    def _got_cell_changed(self):
        if self['got cell'] is True:
            self.region.setBrush((0, 255, 0, 100))
        else:
            self.region.setBrush((0, 0, 255, 100))
        self.region.update()  # pg bug
                