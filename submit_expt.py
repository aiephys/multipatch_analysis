"""
Script used to submit completed experiment to database.
"""

import acq4.util.Canvas, acq4.util.DataManager
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore


class ExperimentSubmitUi(QtGui.QWidget):
    def __init__(self):
        self.path = None
        
        
        QtGui.QWidget.__init__(self)
        
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        
        self.hsplit = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.hsplit, 0, 0)

        self.file_tree = QtGui.QTreeWidget()
        self.hsplit.addWidget(self.file_tree)
        
        self.canvas = acq4.util.Canvas.Canvas()
        self.hsplit.addWidget(self.canvas)
        
    def set_path(self, path):
        self.path = path
        self._reload_file_tree()
        
    def _reload_file_tree(self):
        self.file_tree.clear()
        
        dh = acq4.util.DataManager.getDirHandle(self.path)
        root = self.file_tree.invisibleRootItem()
        self._fill_file_tree(dh, root)
        
    def _fill_file_tree(self, dh, root):
        for fname in dh.ls():
            item = QtGui.QTreeWidgetItem([fname])
            root.addChild(item)
            fh = dh[fname]
            item.fh = fh
            if fh.isDir():
                self._fill_file_tree(fh, item)
        
        
if __name__ == '__main__':
    import sys
    app = pg.mkQApp()
    
    path = sys.argv[1]
    ui = ExperimentSubmitUi()
    ui.resize(1000, 800)
    ui.show()
    ui.set_path(path)
    
    if sys.flags.interactive == 0:
        app.exec_()
