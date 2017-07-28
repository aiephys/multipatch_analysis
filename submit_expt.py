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
        self.file_tree.setSelectionMode(self.file_tree.ExtendedSelection)
        self.hsplit.addWidget(self.file_tree)
        
        self.canvas = acq4.util.Canvas.Canvas()
        self.hsplit.addWidget(self.canvas)
        
        self.hsplit.setSizes([400, 1000])
        
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
            
            if hasattr(item, 'type_selected'):
                item.type_selected.connect(self._item_type_selected)
            
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

    def _item_type_selected(self, item, typ):
        for item in self.file_tree.selectedItems():
            item.set_type(typ)


class ImageTreeItem(pg.TreeWidgetItem):
    
    class Signals(QtCore.QObject):
        type_selected = QtCore.Signal(object, object)
    
    def __init__(self, fh):
        self._sigprox = ImageTreeItem.Signals()
        self.type_selected = self._sigprox.type_selected

        pg.TreeWidgetItem.__init__(self, [fh.shortName(), 'ignore', fh.info()['objective']])
        self.fh = fh
        #self.combo = QtGui.QComboBox()
        self.types = ['ignore', 'type 1', 'type 2']
        #for typ in self.types:
            #self.combo.addItem(typ)
        #self.setWidget(1, self.combo)
        self.menu = QtGui.QMenu()
        for typ in self.types:
            act = self.menu.addAction(typ, self._type_selected)
            
    def _type_selected(self):
        action = self.treeWidget().sender()
        text = str(action.text()).strip()
        self.set_type(text)
        self.type_selected.emit(self, text)
            
    def itemClicked(self, col):
        if col != 1:
            return
        tw = self.treeWidget()
        x = tw.header().sectionPosition(col)
        y = tw.header().height() + tw.visualItemRect(self).bottom()
        self.menu.popup(tw.mapToGlobal(QtCore.QPoint(x, y)))
        return None
    
    def set_type(self, typ):
        self.setText(1, typ)
        




        
if __name__ == '__main__':
    import sys
    app = pg.mkQApp()
    
    path = sys.argv[1]
    ui = ExperimentSubmitUi()
    ui.resize(1300, 800)
    ui.show()
    ui.set_path(path)
    
    if sys.flags.interactive == 0:
        app.exec_()
