from collections import OrderedDict
import numpy as np
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import pyqtgraph.dockarea


class NDSlicer(QtGui.QWidget):
    def __init__(self, axes):
        QtGui.QWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        self.dockarea = pg.dockarea.DockArea()
        self.layout.addWidget(self.dockarea)

        self.viewers = []
        
        self.axes = axes.copy()
        for ax in axes:
            axes[ax].setdefault('index', 0)
        
        self.params = pg.parametertree.Parameter(name='params', type='group', children=[
            {'name': 'index', 'type': 'group', 'children': [{'name': ax, 'type': 'int', 'value': axes[ax]['index']} for ax in axes]},
            {'name': '1D views', 'type': 'group'},
            MultiAxisParam(ndim=2, slicer=self),
        ])

        for ax in axes:
            ch = pg.parametertree.Parameter.create(name=ax, type='bool', value=True)
            self.params.child('1D views').addChild(ch)
            ch.sigValueChanged.connect(self.one_d_show_changed)
            ch.viewer = self.add_view(axes=[ax])
        
        self.ptree = pg.parametertree.ParameterTree(showHeader=False)
        self.ptree.setParameters(self.params, showTop=False)
        self.ptree_dock = pg.dockarea.Dock("view selection")
        self.ptree_dock.addWidget(self.ptree)
        self.dockarea.addDock(self.ptree_dock, 'left')
        
    def add_view(self, axes):
        dock = pg.dockarea.Dock("viewer", area=self.dockarea)
        if len(axes) == 1:
            viewer = OneDViewer(axes)
        elif len(axes) == 2:
            viewer = TwoDViewer(axes)
        dock.addWidget(viewer)
        self.dockarea.addDock(dock, 'right')
        viewer.dock = dock
        viewer.selection_changed.connect(self.viewer_selection_changed)
        self.viewers.append(viewer)
        return viewer

    def one_d_show_changed(self, param):
        param.viewer.dock.setVisible(param.value())
        
    def viewer_selection_changed(self, viewer, axes):
        for ax,val in axes.items():
            self.axes[ax]['index'] = val
            self.params['index', ax] = val
        for viewer in self.viewers:
            viewer.set_index(self.axes)


class MultiAxisParam(pg.parametertree.types.GroupParameter):
    def __init__(self, ndim, slicer):
        self.ndim = ndim
        self.slicer = slicer
        pg.parametertree.types.GroupParameter.__init__(self, name="%dD views"%ndim, addText="Add new..")
        
    def addNew(self):
        axis_names = list(self.slicer.axes.keys())
        param = pg.parametertree.Parameter(name="%dD view" % self.ndim, autoIncrementName=True, type='group', children=[
            {'name': 'axis %d' % i, 'type': 'list', 'values': axis_names, 'value': axis_names[i]} for i in range(self.ndim)
        ], removable=True)
        self.addChild(param)
        viewer = self.slicer.add_view(axes=list(self.slicer.axes.keys())[:self.ndim])
        param.viewer = viewer
        param.sigTreeStateChanged.connect(self.axes_changed)

    def axes_changed(self, param, changes):
        axes = [ch.value() for ch in param.children()]
        param.viewer.set_axes(axes)


class OneDViewer(pg.PlotWidget):
    selection_changed = QtCore.Signal(object, object)  # self, {axis: selection}
    
    def __init__(self, axes):
        pg.PlotWidget.__init__(self)
        self.line = self.addLine(x=0, movable=True)
        self.set_axes(axes)
        self.line.sigDragged.connect(self.line_moved)
        
    def set_axes(self, axes):
        assert len(axes) == 1
        self.axes = axes
        self.setLabels(bottom=axes[0])

    def set_index(self, axes):
        self.line.setValue(axes[self.axes[0]]['index'])

    def line_moved(self):
        self.selection_changed.emit(self, {self.axes[0]: self.line.value()})
        

class TwoDViewer(pg.ImageView):
    selection_changed = QtCore.Signal(object, object)  # self, {axis: selection, ...}

    def __init__(self, axes):
        self.plot = pg.PlotItem()
        pg.ImageView.__init__(self, view=self.plot)
        self.plot.invertY(False)
        self.plot.setAspectLocked(False)
        self.lines = [self.plot.addLine(x=0), self.plot.addLine(y=0)]
        self.set_axes(axes)
        for line in self.lines:
            line.sigDragged.connect(self.line_moved)
        
    def set_axes(self, axes):
        assert len(axes) == 2
        self.axes = axes
        self.plot.setLabels(left=axes[1], bottom=axes[0])
        
    def set_index(self, axes):
        for i,line in enumerate(self.lines):
            line.setValue(axes[self.axes[i]]['index'])

    def line_moved(self):
        axes = {self.axes[i]: self.lines[i].value() for i in (0, 1)}
        self.selection_changed.emit(self, axes)


if __name__ == '__main__':
    pg.mkQApp()
    pg.dbg()
    
    axes = OrderedDict([
        ('X', {}),
        ('Y', {}),
        ('Z', {}),
    ])
        
    nds = NDSlicer(axes)
    nds.show()
                