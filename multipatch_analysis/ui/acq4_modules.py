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

        self.create_json_btn = QtGui.QPushButton("Create Cell Location Json")
        self.layout.addWidget(self.create_json_btn, 1, 0)
        self.create_json_btn.setEnabled(False)
        self.create_json_btn.clicked.connect(self.create_lims_json)

        self.save_json_btn = QtGui.QPushButton("Save and Upload Json")
        self.layout.addWidget(self.save_json_btn, 2, 0)
        self.save_json_btn.setEnabled(False)
        self.save_json_btn.clicked.connect(self.save_json_and_trigger)

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
        self.slice_name = self.base_dir.info()['specimen_ID'].strip()
        self.slice_id = lims.specimen_id_from_name(self.slice_name)
        clusters = lims.cell_cluster_ids(self.slice_id)
        try:
            aff_image_path = lims.specimen_20x_image(self.slice_name)
        except KeyError:
            raise Exception('No Slice Selected')

        if os.path.exists(safe_aff_image_path) == False:
            raise Exception('Couldn\'t find image "%s". Check you have the selected the correct slice.' % safe_aff_image_path)
        
        aff_image_name = os.path.split(aff_image_path)[-1]
        
        save_path = os.path.join(self.base_path, aff_image_name)
        safe_save_path = lims.lims.safe_system_path(save_path)

        if os.path.exists(safe_save_path) == False:
            with pg.BusyCursor():
                shutil.copy2(safe_aff_image_path,safe_save_path)

                self.base_dir.indexFile(aff_image_name)
                print("20x image moved")
        self.image_20 = safe_save_path
        self.create_json_btn.setEnabled(True)
        self.save_json_btn.setEnabled(False)

    def create_lims_json(self):
        """
        Checks that markers are loaded into the canvas
        Gets timestamp
        Creates dictionary of cell locations and enables the save button
        """
        items = self.mosaic_editor.canvas.items 
        markers = [i for i in items if isinstance(i, MarkersCanvasItem)]
        if len(markers) != 1:
            raise Exception("Must have exactly 1 Markers item in the canvas.")
        
        for i in items:
                try:
                    if i.opts['handle'].name() == self.image_20:
                        image = i
                except AttributeError:
                    pass

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
        ts = self.site_dh.info()['__timestamp__']

        
        self.cluster_id = []
        clusters = lims.cell_cluster_ids(self.slice_id)
        for cid in clusters:
            if lims.specimen_metadata(cid) != None:
                cluster_tmsp = lims.specimen_metadata(cid).get('acq_timestamp')
                if cluster_tmsp == ts:
                    self.cluster_id.append(cid)
                    self.cluster_name = lims.specimen_name(cid)
            

        if self.cluster_id == None:
            raise Exception("Couldn't match Multipatch Log to data in LIMS")
        
        if len(self.cluster_id) > 1:
            raise Exception("Multiple clusters found with same timestamp")

        pipette_path = os.path.join(self.site_dh.name(), 'pipettes.yml')
        print(pipette_path)

        if os.path.exists(pipette_path) == False:
            raise Exception("Could not find pipette log")

        pipette_log = yaml.load(open(pipette_path))

        pipettes = markers[0].saveState()['markers']

        self.data = {}  
        self.data['cells'] = []
        for cell in pipettes:
            cell_name = int(cell[0].split('_')[-1])
            p = cell[1]
            x_float = markers[0].graphicsItem().mapToItem(image.graphicsItem(),pg.Point(*p[:2])).x()
            y_float = markers[0].graphicsItem().mapToItem(image.graphicsItem(),pg.Point(*p[:2])).y()
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
        print('Found {} cells in the cluster'.format(len(self.data['cells'])))

        self.save_json_btn.setEnabled(True)

        
    def save_json_and_trigger(self):
        """
        creates the json file and trigger file and moves them to the proper
        incoming folders for LIMS
        disables create json button
        """

        json_save_file = self.cluster_name + "_ephys_cell_cluster.json"

        incoming_path = lims.get_incoming_dir(self.slice_name)
        if incoming_path == None:
            raise Exception("Could not find incoming directory in LIMS")

        incoming_save_path = os.path.sep*2 + os.path.join(*incoming_path.split("/"))
        incoming_json_save_path = os.path.join(incoming_save_path, json_save_file)

        with open(incoming_json_save_path, 'w') as outfile:  
            json.dump(self.data, outfile)


        trigger_file = self.cluster_name + "_ephys_cell_cluster.mpc"

        trigger_path = lims.get_trigger_dir(self.slice_name)
        if trigger_path == None:
            raise Exception("Could not find trigger directory in LIMS")

        trigger_path = os.path.sep*2 + os.path.join(*trigger_path.split("/"))
        
        trigger_save_path = os.path.join(trigger_path,trigger_file)
        with open(trigger_save_path, 'w') as the_file:
            the_file.write("specimen_id: {}\n".format(self.cluster_id[0]))     
            the_file.write("cells_info: '{}'\n".format(incoming_json_save_path))
       
        raise Exception("Thanks for tagging your cells!")
        self.save_json_btn.setEnabled(False)
        

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
