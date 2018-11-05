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
        self.data = None
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

    def set_data(self, data, axes=None):
        self.data = data
        axes = axes or  {}
        for ax in axes:
            self.axes[ax].update(axes[ax])
        for ax in self.axes:
            self.axes[ax]['index'] = self.axes[ax]['values'][0]
        for viewer in self.viewers:
            viewer.set_data(self.data, self.axes)
            viewer.set_index(self.index())
        
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
        viewer.set_data(self.data, self.axes)
        viewer.set_index(self.index())    
        return viewer

    def one_d_show_changed(self, param):
        param.viewer.dock.setVisible(param.value())
        
    def viewer_selection_changed(self, viewer, axes):
        for ax,val in axes.items():
            self.axes[ax]['index'] = val
            self.params['index', ax] = val
        index = self.index()
        for viewer in self.viewers:
            viewer.set_index(index)

    def index(self):
        index = {ax:val['index'] for ax,val in self.axes.items()}
        return index        


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
        param.viewer.set_selected_axes(axes)


class Viewer(object):
    def __init__(self, ax):
        self.data = None
        self.data_axes = None
        self.index = None
        self.set_selected_axes(ax)

    def set_data(self, data, axes):
        self.data = data
        self.data_axes = axes
        self.index = {ax:0 for ax in axes}
        self.update_display()
        
    def set_selected_axes(self, axes):
        self.selected_axes = axes
        if self.index is not None:
            self.set_index(self.index)

    def set_index(self, index):
        self.index = index
        self.update_display()
        
    def update_display(self):
        raise NotImplementedError()
        
    def get_data(self):
        sl = []
        for ax in self.data_axes:
            x = self.index[ax]
            i = self._closest_axis_index(ax, x)
            sl.append(i)
        order = []
        for ax in self.selected_axes:
            i = list(self.data_axes.keys()).index(ax)
            order.append(i)
            sl[i] = slice(None)
        data = self.data[tuple(sl)].transpose(np.argsort(order))
        
        return data

    def _closest_axis_index(self, axis, x):
        vals = self.data_axes[axis]['values']
        return np.argmin(np.abs(vals - x))


class OneDViewer(Viewer, pg.PlotWidget):
    selection_changed = QtCore.Signal(object, object)  # self, {axis: value, ...}
    
    def __init__(self, axes):
        pg.PlotWidget.__init__(self)
        self.line = self.addLine(x=0, movable=True)
        self.curve = self.plot()

        Viewer.__init__(self, axes)

        self.line.sigDragged.connect(self.line_moved)

    def set_index(self, index):
        axis = self.selected_axes[0]
        self.line.setValue(self.data_axes[axis]['index'])
        Viewer.set_index(self, index)

    def line_moved(self):
        self.selection_changed.emit(self, {self.selected_axes[0]: self.line.value()})
        
    def update_display(self):
        if self.data is None:
            self.curve.setData([])
            return
        axis = self.selected_axes[0]
        self.setLabels(bottom=axis)
        data = self.get_data()
        axvals = self.data_axes[axis]['values']
        self.curve.setData(axvals, data)
        

class TwoDViewer(Viewer, pg.ImageView):
    selection_changed = QtCore.Signal(object, object)  # self, {axis: value, ...}

    def __init__(self, axes):
        self.plot = pg.PlotItem()
        
        pg.ImageView.__init__(self, view=self.plot)
        self.plot.invertY(False)
        self.plot.setAspectLocked(False)
        self.lines = [self.plot.addLine(x=0, movable=True), self.plot.addLine(y=0, movable=True)]
        Viewer.__init__(self, axes)
        for line in self.lines:
            line.sigDragged.connect(self.line_moved)
        
    def set_index(self, index):
        for i,line in enumerate(self.lines):
            ax = self.selected_axes[i]
            line.setValue(index[ax])
        Viewer.set_index(self, index)

    def line_moved(self):
        axes = {self.selected_axes[i]: self.lines[i].value() for i in (0, 1)}
        self.selection_changed.emit(self, axes)

    def update_display(self):
        if self.data is None:
            self.setImage(np.zeros((1, 1)))
            return
        axes = self.selected_axes
        self.plot.setLabels(left=axes[1], bottom=axes[0])
        data = self.get_data()
        xvals = self.data_axes[axes[0]]['values']
        yvals = self.data_axes[axes[1]]['values']
        
        scale = [xvals[1]-xvals[0], yvals[1]-yvals[0]]
        self.setImage(data, pos=[xvals[0]-scale[0]*0.5, yvals[0]-scale[1]*0.5], scale=scale)
        


if __name__ == '__main__':
    pg.mkQApp()
    pg.dbg()
    
    data = np.random.normal(size=(10, 20, 40)) * 0.1
    data += np.sin(np.linspace(0, 2*np.pi, data.shape[0]))[:, None, None]
    data += np.sin(np.linspace(0, 4*np.pi, data.shape[1]))[None, :, None]
    data += np.sin(np.linspace(0, 8*np.pi, data.shape[2]))[None, None, :]
    
    axes = OrderedDict([
        ('X', {'values': np.arange(data.shape[0])}),
        ('Y', {'values': np.arange(data.shape[1]) * 200}),
        ('Z', {'values': np.arange(data.shape[2]) * 0.01 + 3.4}),
    ])
        
    nds = NDSlicer(axes)
    nds.show()
    
    nds.set_data(data)
    