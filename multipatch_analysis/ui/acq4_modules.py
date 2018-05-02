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

        self.create_json_btn = QtGui.QPushButton("Save Cell Location Json")
        self.layout.addWidget(self.create_json_btn, 1, 0)
        self.create_json_btn.clicked.connect(self.create_lims_json)


    def load_clicked(self):
        """
        Checks that directory is pointed at a slice
        Check LIMS for 20x biocytin image and if 20x image is saved locally.
        If image not saved locally then copies file from LIMS.
        Opens button to create a JSON of the cell locations
        """
        self.base_dir = self.mosaic_editor.ui.fileLoader.baseDir()
        if self.base_dir.info()['dirType'] != 'Slice':
            raise Exception('No Slice Selected')
        self.base_path = self.base_dir.name()
        self.slice_name = self.base_dir.info()['specimen_ID']
        self.slice_id = lims.specimen_id_from_name(self.slice_name)
        clusters = lims.cell_cluster_ids(self.slice_id)
        """if len(clusters) == 0:
                raise Exception("No Clusters Recorded in LIMS")"""
        try:
            aff_image_path = lims.specimen_20x_image(self.slice_name)
        except KeyError:
            raise Exception('No Slice Selected')
        
        #check the image file path
        safe_aff_image_path = lims.lims.safe_system_path(aff_image_path)

        if os.path.exists(safe_aff_image_path) == False:
            raise Exception("Couldn't find image in LIMS. Check you have the selected the correct slice.")
        
        aff_image_name = os.path.split(aff_image_path)[-1]
        
        jpg_image_name = os.path.splitext(aff_image_name)[0] + ".jpg"

        image_service_url = "http://lims2/cgi-bin/imageservice?path="
        zoom = 6        #need to check that this ACQ4 can support this size

        full_url = image_service_url + aff_image_path + "&zoom={}".format(str(zoom))
        save_path = os.path.join(self.base_path, jpg_image_name)
        safe_save_path = lims.lims.safe_system_path(save_path)

        if os.path.exists(safe_save_path) == False:
            image = urllib.URLopener()
            image.retrieve(full_url, safe_save_path)
            self.base_dir.indexFile(jpg_image_name)
            print("20x image moved")
        self.image_20 = safe_save_path
                

    def create_lims_json(self):
        """
        Checks that markers are loaded into the canvas
        Gets timestamp
        Creates json for each cell clusters and saves locally
        copies json to correct incoming folder
        creates trigger file and drops in correct incoming folder
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

        
        cluster_id = []
        clusters = lims.cell_cluster_ids(self.slice_id)
        for cid in clusters:
            cluster_tmsp = lims.specimen_metadata(cid)['acq_timestamp']
            if cluster_tmsp == ts:
                cluster_id.append(cid)
                cluster_name = lims.specimen_name(cid)
            

        if cluster_id == None:
            raise Exception("Couldn't match Multipatch Log to data in LIMS")
        
        if len(cluster_id) > 1:
            raise Exception("Multiple clusters found with same timestamp")

        pipette_path = os.path.join(self.site_dh.name(), 'pipettes.yml')

        if os.path.exists(pipette_path) == False:
            raise Exception("Could not find pipette log")

        pipette_log = yaml.load(open(pipette_path))

        pipettes = markers[0].saveState()['markers']

        data = {}  
        data['cells'] = []
        for cell in pipettes:
            cell_name = int(cell[0].split('_')[-1])
            p = cell[1]
            data['cells'].append({      
                'ephys_cell_id': cell_name,                  
                'coordinates_20x': 
                {
                "x": markers[0].graphicsItem().mapToItem(image.graphicsItem(),pg.Point(*p[:2])).x(), 
                "y": markers[0].graphicsItem().mapToItem(image.graphicsItem(),pg.Point(*p[:2])).y()
                },      
                'ephys_qc_result': ('pass' if pipette_log[cell_name]['got_data'] == 'true' else 'fail'),                 
                'start_time_sec': time.mktime(pipette_log[cell_name]['patch_start'].timetuple())
            })
        
        # as defined by requirements from technology

        json_save_file = cluster_name + "_ephys_cell_cluster.json"

        incoming_path = lims.get_incoming_dir(self.slice_name)
        if incoming_path == None:
            raise Exception("Could not find incoming directory in LIMS")

        incoming_path = os.path.sep*2 + os.path.join(*incoming_path.split("/"))
        local_json_save_path = os.path.join(parent.name(), json_save_file)

        with open(local_json_save_path, 'w') as outfile:  
            json.dump(data, outfile)

        incoming_json_path = os.path.join(incoming_path, json_save_file)
        
        # uncomment line below to copy json to incoming folder
        #shutil.copy2(local_json_save_path, incoming_json_path)

        trigger_file = cluster_name + "_ephys_cell_cluster.mpc"

        trigger_path = lims.get_trigger_dir(self.slice_name)
        if trigger_path == None:
            raise Exception("Could not find trigger directory in LIMS")

        trigger_path = os.path.sep*2 + os.path.join(*trigger_path.split("/"))
        
        # prevents dropping trigger file in incoming folder
        # comment out line below to send trigger file to LIMS
        trigger_path = os.path.join(parent.name(),trigger_file)
        with open(trigger_path, 'w') as the_file:
            the_file.write("specimen_id: {}\n".format(cluster_id[0]))     #verify if this is the cluster specimen or the slice specimen
            the_file.write("cells_info: '{}'\n".format(incoming_json_path))
        raise Exception('Trigger File Saved Locally')
        


MosaicEditor.addExtension("Multi Patch", {
    'type': 'ctrl',
    'builder': MultiPatchMosaicEditorExtension,
    'pos': ('top', 'Canvas'),
    'size': (600, 200),
})
