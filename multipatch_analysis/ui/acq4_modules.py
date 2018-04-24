import acq4.pyqtgraph as pg
from acq4.pyqtgraph.Qt import QtCore, QtGui
from acq4.modules.Module import Module
from acq4.Manager import getManager
# from acq4.analysis.modules.MosaicEditor import MosaicEditor
# from acq4.util.Canvas.items.CanvasItem import CanvasItem
from acq4.modules.MosaicEditor import MosaicEditor
from acq4.util.Canvas.items.MarkersCanvasItem import MarkersCanvasItem 
from acq4.util.Canvas.items.MultiPatchLogCanvasItem import MultiPatchLogCanvasItem 
from . import submit_expt
from . import multipatch_nwb_viewer
from .. import lims
import os
from acq4.util import FileLoader
import shutil
import json
import urllib

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

        self.load_btn = QtGui.QPushButton("Load 20x")
        self.layout.addWidget(self.load_btn, 0, 0)
        self.load_btn.clicked.connect(self.load_clicked)

        self.create_json_btn = QtGui.QPushButton("Save Cell Location Json")
        self.layout.addWidget(self.create_json_btn, 1, 0)
        self.create_json_btn.clicked.connect(self.create_lims_json)

        self.create_test_btn = QtGui.QPushButton('Testing Button')
        self.layout.addWidget(self.create_test_btn, 2, 0)
        self.create_test_btn.clicked.connect(self.test_button)

    def test_button(self):
        """place to test things
        """
        items = self.mosaic_editor.canvas.items 
        logs = [i for i in items if isinstance(i, MultiPatchLogCanvasItem)]
        print logs[0].saveState()['currentTime']
        

    def load_clicked(self):
        """
        Checks that directory is pointed at a slice
        Check LIMS for 20x biocytin image and if 20x image is saved locally.
        If image not saved locally then copies file from LIMS.
        Opens button to create a JSON of the cell locations
        """
        self.base_dir = self.mosaic_editor.ui.fileLoader.baseDir()
        if self.base_dir.info()['dirType'] == 'Slice':
            self.base_path = self.base_dir.name()
            self.slice_name = self.base_dir.info()['specimen_ID']
            self.slice_id = lims.specimen_id_from_name(self.slice_name)
            clusters = lims.cell_cluster_ids(self.slice_id)
            """if len(clusters) == 0:
                                                    raise Exception("No Clusters Recorded in LIMS")"""
            try:
                aff_image_path = lims.specimen_20x_image(self.slice_name)
                #check the image file path
                safe_aff_image_path = lims.lims.safe_system_path(aff_image_path)

                if os.path.exists(safe_aff_image_path) == True:
                    print('File Exists in LIMS')
                    aff_image_name = aff_image_path.split("/")[-1]
                    jpg_image_name = aff_image_name.split(".aff")[0] + ".jpg"

                    image_service_url = "http://lims2/cgi-bin/imageservice?path="
                    zoom = 5        #need to check that this ACQ4 can support this size

                    full_url = image_service_url + aff_image_path + "&zoom={}".format(str(zoom))
                    print(full_url)
                    save_path = self.base_path + '/' + jpg_image_name
                    safe_save_path = lims.lims.safe_system_path(save_path)

                    if os.path.exists(safe_save_path) == True:
                        print("20x already moved")

                    elif os.path.exists(safe_save_path) == False:
                        image = urllib.URLopener()
                        image.retrieve(full_url, safe_save_path)
                        print("20x image moved")
                                             
                else:
                    raise Exception("Couldn't find image in LIMS. Check you have the selected the correct slice.")
                
            except KeyError:
                raise Exception('No Slice Selected')
        else:
            raise Exception('No Slice Selected')
       


    def create_lims_json(self):
        """
        Creates json for each cell clusters and saves locally
        copies json to correct incoming folder
        creates trigger file and drops in correct incoming folder
        """
        
        items = self.mosaic_editor.canvas.items 

        markers = [i for i in items if isinstance(i, MarkersCanvasItem)]
        if len(markers) != 1:
            raise Exception("Must have exactly 1 Markers item in the canvas.")
        markers = markers[0].saveState()['markers']

        logs = [i for i in items if isinstance(i, MultiPatchLogCanvasItem)]
        if len(logs) != 1:
            raise Exception("Must have exactly 1 Log item in the canvas.")
        multipatch_tmsp = logs[0].saveState()['currentTime']

        cluster_name = '01'
        cluster_id = []
        clusters = lims.cell_cluster_ids(self.slice_id)
        for cid in clusters:
            cluster_tmsp = lims.specimen_metadata[cid]['acq_timestamp']
            if cluster_tmsp == multipatch_tmsp:
                cluster_id = cid
                cluster_name = lims.specimen_name(cid).split(".")[-1]
                break

        """if cluster_id == []:
                                    raise Exception("Couldn't match Multipatch Log to data in LIMS")"""
            
        data = {}  
        data['cells'] = []
        for cell in markers:
            data['cells'].append({      
                'ephys_cell_id': cell[0],                        
                'coordinates_20x': {"x": cell[1][0], "y": cell[1][1]},      
                'ephys_qc_result': 'pass',                  #find this value from luke
                'start_time_sec': 345.4                     #find this value from luke
            })
        
        # as defined by requirements from technology

        json_save_file = self.slice_name + "_ephys_cell_cluster_" + cluster_name + ".json"
        print (json_save_file)

        if lims.is_mouse(self.slice_name) == 'mouse':
            incoming_path = "/allen/programs/celltypes/production/incoming/mousecelltypes/"
        
        if lims.is_mouse(self.slice_name) == 'human':
            incoming_path = "/allen/programs/celltypes/production/incoming/humnacelltypes/"
        
        local_json_save_path = os.path.join(self.base_path, json_save_file)

        with open(local_json_save_path, 'w') as outfile:  
            json.dump(data, outfile)

        incoming_json_path = os.path.join(incoming_path, json_save_file)
        #shutil.copy2(local_json_save_path, json_path)
        print (incoming_json_path)

        trigger_file = self.slice_name + "_ephys_cell_cluster_" + cluster_name + "_cells.mpc"
        trigger_path = os.path.join(incoming_path, 'trigger', trigger_file)
        trigger_path = lims.lims.safe_system_path(trigger_path)
        
        
        #prevents dropping trigger file in incoming folder
        trigger_path = os.path.join(self.base_dir.name(),trigger_file)
        print (trigger_path)
        with open(trigger_path, 'w') as the_file:
            the_file.write("specimen_id: {}\n".format(self.slice_id))     #verify if this is the cluster specimen or the slice specimen
            the_file.write("cells_info: '{}'\n".format(incoming_json_path))
        
        


MosaicEditor.addExtension("Multi Patch", {
    'type': 'ctrl',
    'builder': MultiPatchMosaicEditorExtension,
    'pos': ('top', 'Canvas'),
    'size': (600, 200),
})
