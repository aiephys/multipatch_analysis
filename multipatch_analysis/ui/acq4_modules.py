import acq4.pyqtgraph as pg
from acq4.pyqtgraph.Qt import QtCore, QtGui
from acq4.modules.Module import Module
from acq4.Manager import getManager
import acq4.util.DataManager as DataManager

# from acq4.analysis.modules.MosaicEditor import MosaicEditor
# from acq4.util.Canvas.items.CanvasItem import CanvasItem
from acq4.modules.MosaicEditor import MosaicEditor
from acq4.util.Canvas.items.MarkersCanvasItem import MarkersCanvasItem 
from acq4.util.Canvas.items.MultiPatchLogCanvasItem import MultiPatchLogCanvasItem 
from . import submit_expt
from . import multipatch_nwb_viewer
from .. import lims
import os
import shutil
import json
import urllib
from datetime import datetime
import time
import yaml
import affpyramid

from acq4.util.Canvas.items.CanvasItem import CanvasItem
from acq4.util.Canvas.items import registerItemType
from . import dashboard


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

        self.load_btn = QtGui.QPushButton("Download 20x")
        self.layout.addWidget(self.load_btn, 0, 0)
        self.load_btn.clicked.connect(self.load_clicked)

        self.submit_btn = pg.FeedbackButton("Submit cell positions to LIMS")
        self.layout.addWidget(self.submit_btn, 1, 0)
        self.submit_btn.setEnabled(False)
        self.submit_btn.clicked.connect(self.submit)

    def load_clicked(self):
        """
        Checks that directory is pointed at a slice
        Check LIMS for 20x biocytin image and if 20x image is saved locally.
        If image not saved locally then copies file from LIMS.
        Enables button to create a JSON of the cell locations
        """
        self.base_dir = self.mosaic_editor.ui.fileLoader.baseDir()
        if self.base_dir.info()['dirType'] != 'Slice':
            raise Exception('No Slice Selected')
        self.base_path = self.base_dir.name()
        self.specimen_name = self.base_dir.info()['specimen_ID'].strip()

        aff_image_path = lims.specimen_20x_image(self.specimen_name, treatment='Biocytin')
        if aff_image_path is None:
            raise Exception("No 20x biocytin image for specimen %s" % self.specimen_name)

        if not os.path.exists(aff_image_path):
            raise Exception("No image found at path '%s'" % aff_image_path)
        
        aff_image_name = os.path.split(aff_image_path)[-1]
        save_path = os.path.join(self.base_path, aff_image_name)

        if not os.path.exists(save_path):
            with pg.BusyCursor():
                shutil.copy2(aff_image_path, save_path)
                self.base_dir.indexFile(aff_image_name)
        
        self.image_20 = save_path
        self.submit_btn.setEnabled(True)

    def submit(self):
        try:
            self.create_lims_json()
            self.save_json_and_trigger()
            self.submit_btn.success()
        except:
            self.submit_btn.failure()
            raise

    def create_lims_json(self):
        """
        Checks that markers are loaded into the canvas
        Gets timestamp
        Creates dictionary of cell locations and enables the save button
        """
        items = self.mosaic_editor.canvas.items 

        # Find item that marks cell positions
        markers = [i for i in items if isinstance(i, MarkersCanvasItem)]
        if len(markers) != 1:
            raise Exception("Must have exactly 1 Markers item in the canvas.")
        
        # Find the 20x image that was loaded in to the editor
        image = None
        for i in items:
            try:
                if i.opts['handle'].name() == self.image_20:
                    image = i
            except AttributeError:
                pass
        if image is None:
            raise Exception("Could not find 20x image loaded into mosaic editor.")

        # Infer the site folder to submit from based on the files that have been loaded.
        # This is really indirect, but it seems to work.
        sites = set()
        for item in items:
            handle = item.opts.get('handle')
            if handle is None:
                continue
            parent = handle.parent()
            if parent.shortName().startswith('site_'):
                sites.add(parent)
        if len(sites) == 0:
            raise Exception("No files loaded from site; can't determine which site to submit.")
        if len(sites) > 1:
            raise Exception("Files loaded from multiple sites; can't determine which site to submit.")
        self.site_dh = list(sites)[0]

        # Grab the timestamp for this site; we need it to look up the cell cluster from LIMS
        ts = self.site_dh.info()['__timestamp__']

        # Find the LIMS CellCluster ID for this experiment 
        cluster_ids = lims.expt_cluster_ids(self.specimen_name, ts)
        if len(cluster_ids) == 0:
            raise Exception("No CellCluster found for %s %s" % (self.specimen_name, ts))
        if len(cluster_ids) > 1:
            raise Exception("Multiple CellClusters found for %s %s" % (self.specimen_name, ts))
        self.cluster_id = cluster_ids[0]
        self.cluster_name = lims.specimen_name(self.cluster_id)

        # Load up the piettes.yml file so we know which cells are considered "ephys pass"
        pipette_path = os.path.join(self.site_dh.name(), 'pipettes.yml')
        if not os.path.exists(pipette_path):
            raise Exception("Could not find pipettes.yml")
        pipette_log = yaml.load(open(pipette_path))

        # grab info about cell positions from Markers
        pipettes = markers[0].saveState()['markers']

        # Construct JSON to send to LIMS.
        self.data = {}  
        self.data['cells'] = []
        for cell in pipettes:
            cell_name = int(cell[0].split('_')[-1])
            p = cell[1]
            x_float = markers[0].graphicsItem().mapToItem(image.graphicsItem(), pg.Point(*p[:2])).x()
            y_float = markers[0].graphicsItem().mapToItem(image.graphicsItem(), pg.Point(*p[:2])).y()
            self.data['cells'].append({      
                'ephys_cell_id': cell_name,                  
                'coordinates_20x': 
                {
                "x": int(round(x_float)),   # round the float and change to integer for technology requirements
                "y": int(round(y_float))    # round the float and change to integer for technology requirements
                },      
                'ephys_qc_result': ('pass' if pipette_log[cell_name]['got_data'] == True else 'fail'),                 
                'start_time_sec': time.mktime(pipette_log[cell_name]['patch_start'].timetuple())
            })
        
    def save_json_and_trigger(self):
        """
        creates the json file and trigger file and moves them to the proper
        incoming folders for LIMS
        disables create json button
        """

        json_save_file = 'synphys_' + self.cluster_name + "_ephys_cell_cluster.json"

        # Ask LIMS where to put incoming files
        incoming_path = lims.get_incoming_dir(self.specimen_name)
        if incoming_path is None:
            raise Exception("Could not find LIMS incoming directory for specimen '%s'" % self.specimen_name)

        incoming_save_path = '//' + incoming_path.lstrip('/')
        if not os.path.isdir(incoming_save_path):
            raise Exception("LIMS incoming directory for specimen '%s' does not exist: '%s'" % (self.specimen_name, incoming_save_path))
        
        # Save payload json file
        incoming_json_save_path = incoming_save_path + '/' + json_save_file # using '/' on windows because it is a network path
        with open(incoming_json_save_path, 'w') as outfile:  
            json.dump(self.data, outfile)

        trigger_file = 'synphys_' + self.cluster_name + "_ephys_cell_cluster.mpc"

        trigger_path = lims.get_trigger_dir(self.specimen_name)
        if trigger_path is None:
            raise Exception("Could not find LIMS trigger directory for specimen '%s'" % self.specimen_name)
        trigger_path = '//' + trigger_path.lstrip('/')
        
        trigger_save_path = os.path.join(trigger_path, trigger_file)
        print trigger_save_path
        print incoming_json_save_path
        with open(trigger_save_path, 'w') as the_file:
            the_file.write("specimen_id: {}\n".format(self.cluster_id))     
            the_file.write("cells_info: '{}'\n".format(incoming_json_save_path))
        

MosaicEditor.addExtension("Multi Patch", {
    'type': 'ctrl',
    'builder': MultiPatchMosaicEditorExtension,
    'pos': ('top', 'Canvas'),
    'size': (600, 200),
})


class DashboardModule(Module):
    """ACQ module for monitoring pipeline status
    """
    moduleDisplayName = "MP Dashboard"
    moduleCategory = "Analysis"

    def __init__(self, manager, name, config):
        Module.__init__(self, manager, name, config)
        self.ui = dashboard.Dashboard(limit=config.get('limit', None), filter_defaults=config.get('filters', None))
        self.ui.resize(1600, 900)
        self.ui.show()
        
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
        opts['defaultUserTransform'] = {'scale': (0.36e-6, 0.36e-6)}            
        CanvasItem.__init__(self, self.affitem, **opts)

    @classmethod
    def checkFile(cls, fh):
        name = fh.shortName()
        if name.endswith('.aff'):
            return 10
        else:
            return 0

registerItemType(AffPyramidCanvasItem)
