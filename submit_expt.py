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

        self.file_tree = pg.TreeWidget()
        self.file_tree.setColumnCount(3)
        self.file_tree.setHeaderLabels(['file', 'category', 'metadata'])
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
            fh = dh[fname]
            item = self._make_item(fh)
            root.addChild(item)
            item.fh = fh
            if fh.isDir():
                self._fill_file_tree(fh, item)
        
    def _make_item(self, fh):
        info = fh.info()
        objtyp = info.get('__object_type__')
        if objtyp in ['ImageFile', 'MetaArray']:
            return ImageTreeItem(fh)
        else:
            item = pg.TreeWidgetItem([fh.shortName()])
            print objtyp, '\t', fh.name()
            return item


class ImageTreeItem(pg.TreeWidgetItem):
    def __init__(self, fh):
        pg.TreeWidgetItem.__init__(self, [fh.shortName(), '', fh.info()['objective']])
        self.fh = fh
        self.combo = QtGui.QComboBox()
        self.types = ['ignore', 'type 1', 'type 2']
        for typ in self.types:
            self.combo.addItem(typ)
        self.setWidget(1, self.combo)
        




        
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
