import acq4.pyqtgraph as pg
from acq4.pyqtgraph.Qt import QtCore, QtGui
from acq4.modules.Module import Module
from acq4.Manager import getManager
# from acq4.analysis.modules.MosaicEditor import MosaicEditor
# from acq4.util.Canvas.items.CanvasItem import CanvasItem
from acq4.modules.MosaicEditor import MosaicEditor
from . import submit_expt
from . import multipatch_nwb_viewer


class MultipatchSubmissionModule(Module):
    """Allows multipatch data submission UI to be invoked as an ACQ4 module.
    
    This is primarily to ensure that metadata changes made via ACQ4 are immediately
    available to (and do noty collide with) the submission tool, and vice-versa. 
    """
    moduleDisplayName = "MP Submission Tool"
    moduleCategory = "Analysis"

    def __init__(self, manager, name, config):
        Module.__init__(self, manager, name, config)
        self.ui = submit_expt.ExperimentSubmitUi()
        self.ui.resize(1600, 900)
        self.ui.show()
        
        self.load_from_dm_btn = QtGui.QPushButton("load from data manager")
        self.ui.left_layout.insertWidget(0, self.load_from_dm_btn)
        self.load_from_dm_btn.clicked.connect(self.load_from_dm_clicked)
        
    def load_from_dm_clicked(self):
        man = getManager()
        sel_dir = man.currentFile
        self.ui.set_path(sel_dir)
        
    def window(self):
        return self.ui

        
class NWBViewerModule(Module):
    """ACQ module for browsing data in NWB files.
    """
    moduleDisplayName = "MP NWB Viewer"
    moduleCategory = "Analysis"

    def __init__(self, manager, name, config):
        Module.__init__(self, manager, name, config)
        self.ui = multipatch_nwb_viewer.MultipatchNwbViewer()
        self.ui.resize(1600, 900)
        self.ui.show()
        
        self.load_from_dm_btn = QtGui.QPushButton("load from data manager")
        self.ui.vsplit.insertWidget(0, self.load_from_dm_btn)
        self.load_from_dm_btn.clicked.connect(self.load_from_dm_clicked)
        
    def load_from_dm_clicked(self):
        man = getManager()
        filename = man.currentFile.name()
        nwb = self.ui.load_nwb(filename)
        
    def window(self):
        return self.ui


class MultiPatchMosaicEditorExtension(QtGui.QWidget):
    def __init__(self, mosaic_editor):
        self.mosaic_editor = mosaic_editor

        QtGui.QWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.label = QtGui.QLabel("OHai.")
        self.layout.addWidget(self.label)

MosaicEditor.addExtension("Multi Patch", {
    'type': 'ctrl',
    'builder': MultiPatchMosaicEditorExtension,
    'pos': ('top', 'Canvas'),
    'size': (600, 200),
})
