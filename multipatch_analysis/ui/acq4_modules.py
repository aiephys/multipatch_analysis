from acq4.pyqtgraph.Qt import QtCore, QtGui
from acq4.modules.Module import Module
from acq4.Manager import getManager
from acq4.util.Canvas.items.CanvasItem import CanvasItem
from acq4.util.Canvas.items import registerItemType
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


class AffPyramidCanvasItem(CanvasItem):
    """For displaying AFF image pyramids
    """
    _typeName = "AFF Image Pyramid"
    
    def __init__(self, handle, **kwds):
        from affpyramid.ui import AffImageItem
        kwds.pop('viewRect', None)
        self.affitem = AffImageItem(handle.name())
        opts = {'movable': True, 'rotatable': True, 'handle': handle}
        opts.update(kwds)
        if opts.get('name') is None:
            opts['name'] = handle.shortName()            
        CanvasItem.__init__(self, self.affitem, **opts)

    @classmethod
    def checkFile(cls, fh):
        name = fh.shortName()
        if name.endswith('.aff'):
            return 10
        else:
            return 0

registerItemType(AffPyramidCanvasItem)
