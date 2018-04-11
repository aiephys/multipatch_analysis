import acq4.pyqtgraph as pg
from acq4.pyqtgraph.Qt import QtCore, QtGui
from acq4.modules.Module import Module
from acq4.Manager import getManager
# from acq4.analysis.modules.MosaicEditor import MosaicEditor
# from acq4.util.Canvas.items.CanvasItem import CanvasItem
from acq4.modules.MosaicEditor import MosaicEditor
from . import submit_expt
from . import multipatch_nwb_viewer
from .. import lims
import os
from acq4.util import FileLoader
import shutil
import json


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

    def load_clicked(self):
        base_dir = self.mosaic_editor.ui.fileLoader.baseDir()
        try:
            spec_id = base_dir.info()['specimen_ID']
            print(spec_id)
            image_path = lims.specimen_20x_image(spec_id)
            print(image_path)
        
            #check the image file path
            if os.path.exists(image_path) == True:
                print('File Exists in LIMS')

                image_name = image_path.split("\\")[-1]
                # need to copy image to current directory
                # need to find current directory
                # hardcode for the time being
                save_path = "/Users/aarono/Scripts/example_data/2017.09.06_000/slice_001"
                save_path = save_path + '/' + image_name
                safe_save_path = lims.lims.safe_system_path(save_path)

                if os.path.exists(safe_save_path) == False:
                    shutil.copy2(image_path, safe_save_path)
                    
                if os.path.exists(safe_save_path) == True:
                    print("20x already moved")

                self.create_json_btn = QtGui.QPushButton("Save Cell Location Json")
                self.layout.addWidget(self.create_json_btn, 1, 0)

                self.create_json_btn.clicked.connect(self.create_lims_json)

                
                
                
            else:
                print("Couldn't find image in LIMS. Check you have the selected the correct slice.")
            #print (base_dir.info())
            #print(os.path.join(base_dir, image_name))
        except KeyError:
            print('No Slice Selected')
        #raise Exception()

    def upload_clicked(self):
        print ("uploaded to LIMS")

        # create trigger file in proper folder
        # name spec_id_ephys_cluster.mpc
        # Chat-IRES-Cre-neo;Ai14-316041.04.01_ephys_cell_cluster_1234_cells.mpc
        # contains specimen_id: <LIMS id of cell cluster>
        # cells_info: '<full path of json file in incoming directory>'
        # specimen_id: 587022731
        # cells_info: '/allen/programs/celltypes/production/incoming/mousecelltypes/Chat-IRES-Cre-neo;Ai14-316041.04.01_ephys_cell_cluster_1234_cells.json'
    
        # either mouse cell types or human cell types
        # /allen/programs/celltypes/production/incoming/mousecelltypes/trigger/ or /allen/programs/celltypes/production/incoming/humancelltypes/trigger/
        
        # need to get cluster id, species, spec_id, spec_name

    def create_lims_json(self):
        # save locally so can be accessed in acq4 for second opinion and then upload function will make a copy
        # that will be dropped in incoming folder

        #stand in until I can move the variables between functions
        base_dir = self.mosaic_editor.ui.fileLoader.baseDir()
        spec_id = base_dir.info()['specimen_ID']
        print(spec_id)

        # need to get cluster id cell locations
        # i don't know what the start time variable is
        data = {}  
        data['cells'] = []  
        data['cells'].append({  
            'ephys_cell_id': 23,
            'coordinates_20x': {"x": 20, "y": 40},
            'ephys_qc_result': 'pass',
            'start_time_sec': 345.4
        })
        data['cells'].append({  
            'ephys_cell_id': 24,
            'coordinates_20x': {"x": 23, "y": 42},
            'ephys_qc_result': 'pass',
            'start_time_sec': 345.4
        })
        data['cells'].append({  
            'ephys_cell_id': 25,
            'coordinates_20x': {"x": 26, "y": 4},
            'ephys_qc_result': 'fail',
            'start_time_sec': 7654.2
        })
        # as defined by requirements from technology
        save_file = spec_id + '_ephys_cell_cluster_1234_cells.json'
        with open(save_file, 'w') as outfile:  
            json.dump(data, outfile)

        self.upload_btn = QtGui.QPushButton("Upload Tags to LIMS")
        self.layout.addWidget(self.upload_btn, 1, 0)

        self.upload_btn.clicked.connect(self.upload_clicked)


MosaicEditor.addExtension("Multi Patch", {
    'type': 'ctrl',
    'builder': MultiPatchMosaicEditorExtension,
    'pos': ('top', 'Canvas'),
    'size': (600, 200),
})
