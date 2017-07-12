from collections import OrderedDict
import os, sys, subprocess
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui


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
