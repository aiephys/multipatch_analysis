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

        self.ctrl_widget = QtGui.QWidget()
        self.ctrl_layout = QtGui.QGridLayout()
        self.ctrl_widget.setLayout(self.ctrl_layout)
        self.hsplit.addWidget(self.ctrl_widget)
        
        self.file_tree = FileTreeWidget()
        self.ctrl_layout.addWidget(self.file_tree, 0, 0, 1, 2)
        
        self.load_btn = QtGui.QPushButton('load files')
        self.load_btn.clicked.connect(self.load_clicked)
        self.ctrl_layout.addWidget(self.load_btn, self.ctrl_layout.rowCount(), 0)
        
        self.canvas = acq4.util.Canvas.Canvas(allowTransforms=False)
        self.hsplit.addWidget(self.canvas)
        
        self.hsplit.setSizes([600, 700])
        
    def set_path(self, path):
        self.path = path
        self.file_tree.set_path(path)

    def load_clicked(self):
        sel = self.file_tree.selectedItems()
        for item in sel:
            fh = item.fh
            self.canvas.addFile(fh)


class FileTreeWidget(pg.TreeWidget):
    def __init__(self):
        pg.TreeWidget.__init__(self)
        self.path = None
        self.setColumnCount(3)
        self.setHeaderLabels(['file', 'category', 'metadata'])
        self.setSelectionMode(self.ExtendedSelection)

    def set_path(self, path):
        self.path = path
        self._reload_file_tree()
        
    def _reload_file_tree(self):
        self.clear()
        
        dh = acq4.util.DataManager.getDirHandle(self.path)
        root = self.invisibleRootItem()
        self._fill_tree(dh, root)
        
    def _fill_tree(self, dh, root):
        for fname in dh.ls():
            fh = dh[fname]
            item = self._make_item(fh)
            item.setExpanded(True)
            if hasattr(item, 'type_selected'):
                item.type_selected.connect(self._item_type_selected)
            
            root.addChild(item)
            item.fh = fh
            if fh.isDir():
                self._fill_tree(fh, item)
        
        for i in range(3):
            self.resizeColumnToContents(i)
        
    def _make_item(self, fh):
        info = fh.info()
        objtyp = info.get('__object_type__')
        if objtyp in ['ImageFile', 'MetaArray']:
            return ImageTreeItem(fh)
        elif fh.shortName().lower().endswith('.nwb'):
            return NwbTreeItem(fh)
        else:
            item = pg.TreeWidgetItem([fh.shortName()])
            print objtyp, '\t', fh.name()
            return item

    def _item_type_selected(self, item, typ):
        for item in self.selectedItems():
            item.set_type(typ)


class TypeSelectItem(pg.TreeWidgetItem):
    """TreeWidgetItem with a type selection menu in the second column.
    """
    class Signals(QtCore.QObject):
        type_selected = QtCore.Signal(object, object)
    
    def __init__(self, fh, types, current_type):
        self.fh = fh
        self._sigprox = ImageTreeItem.Signals()
        self.type_selected = self._sigprox.type_selected
        self.types = types
        pg.TreeWidgetItem.__init__(self, [fh.shortName(), current_type, ''])

        self.menu = QtGui.QMenu()
        for typ in self.types:
            act = self.menu.addAction(typ, self._type_selected)

    def _type_selected(self):
        action = self.treeWidget().sender()
        text = str(action.text()).strip()
        self.set_type(text)
        self.type_selected.emit(self, text)
            
    def set_type(self, typ):
        self.setText(1, typ)
        if typ == 'ignore':
            self.setBachground(1, pg.mkColor(0.8))

    def itemClicked(self, col):
        if col != 1:
            return
        tw = self.treeWidget()
        x = tw.header().sectionPosition(col)
        y = tw.header().height() + tw.visualItemRect(self).bottom()
        self.menu.popup(tw.mapToGlobal(QtCore.QPoint(x, y)))
        return None
        

class NwbTreeItem(TypeSelectItem):
    def __init__(self, fh):
        types = ['ignore', 'MIES physiology']
        TypeSelectItem.__init__(self, fh, types, 'ignore')        
    

class ImageTreeItem(TypeSelectItem):
    def __init__(self, fh):
        info = fh.info()
        meta = info['objective']

        types = ['ignore', 'slice anatomy', 'slice quality stack', 'recording site']
        TypeSelectItem.__init__(self, fh, types, 'ignore')        
        
        self.setText(2, meta)
        colors = info.get('illumination', {}).keys()
        if len(colors) == 0:
            color = 'w'
        elif len(colors) > 1:
            color = 'y'
        else:
            color = {'infrared': (255, 200, 200), 'green': (200, 255, 200), 'blue': (200, 200, 255), 'uv': (240, 200, 255)}[colors[0]]
        self.setBackground(2, pg.mkColor(color))
            




        
if __name__ == '__main__':
    import sys
    app = pg.mkQApp()
    pg.dbg()
    
    path = sys.argv[1]
    ui = ExperimentSubmitUi()
    ui.resize(1300, 800)
    ui.show()
    ui.set_path(path)
    
    if sys.flags.interactive == 0:
        app.exec_()
