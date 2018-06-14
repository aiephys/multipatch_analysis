from collections import OrderedDict
import pyqtgraph as pg


class ExperimentActions(pg.QtCore.QObject):
    """A set of QActions that can be performed on an experiment; used for providing
    consistent context menus across UI tools.
    """
    def __init__(self):
        pg.QtCore.QObject.__init__(self)

        self._experiment = None

        actions = [
            ('Data Manager', self.data_manager),
            ('NWB Viewer', self.nwb_viewer),
            ('Connection Detection', self.connection_detection),
            ('Submission Tool', self.submission_tool),
        ]
        self.actions = OrderedDict()
        for name, callback in actions:
            action = pg.QtGui.QAction(name, self)
            action.triggered.connect(callback)
            self.actions[name] = action

    @property
    def experiment(self):
        """The currently active / selected experiment.

        All triggered actions are performed on this experiment.

        This propery can be overridden to customize lookup of the
        selected experiment.
        """
        return self._experiment

    @experiment.setter
    def experiment(self, expt):
        self._experiment = expt

    def data_manager(self):
        expt = self.experiment
        from acq4.Manager import getManager
        manager = getManager()
        mod = manager.getModule('Data Manager')
        mod.selectFile(expt.path)

    def nwb_viewer(self):
        print(self.experiment)

    def connection_detection(self):
        print(self.experiment)

    def submission_tool(self):
        print(self.experiment)


class CellActions(pg.QtCore.QObject):
    """A set of QActions that can be performed on a cell; used for providing
    consistent context menus across UI tools.
    """
    def __init__(self):
        pg.QtCore.QObject.__init__(self)


class PairActions(pg.QtCore.QObject):
    """A set of QActions that can be performed on a cell pair; used for providing
    consistent context menus across UI tools.
    """
    def __init__(self):
        pg.QtCore.QObject.__init__(self)
