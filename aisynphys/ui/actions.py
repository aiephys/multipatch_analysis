import os, subprocess
from collections import OrderedDict
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import aisynphys.config as config
import aisynphys


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
            ('LIMS Drawing Tool', self.lims_drawing_tool),
            ('Edit pipettes.yml', self.edit_pipettes_yml),
            ('Pair Analysis', self.pair_analysis),
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
        path = os.path.join(os.path.dirname(aisynphys.__file__), '..', 'tools', 'mies_nwb_viewer.py')
        subprocess.Popen('python "%s" "%s"' % (path, self.experiment.nwb_file), shell=True)

    def connection_detection(self):
        path = os.path.join(os.path.dirname(aisynphys.__file__), '..', 'tools', 'connection_detection.py')
        subprocess.Popen('python "%s" "%s"' % (path, self.experiment.nwb_file), shell=True)

    def pair_analysis(self):
        path = os.path.join(os.path.dirname(aisynphys.__file__), '..', 'tools', 'pair_analysis.py')
        subprocess.Popen('python "%s" "%0.3f"' % (path, self.experiment.timestamp), shell=True)

    def submission_tool(self):
        from acq4.Manager import getManager 
        man = getManager()
        st = man.getModule('MultipatchSubmissionModule') 
        st.ui.set_path(man.dirHandle(self.experiment.path))

    def lims_drawing_tool(self):
        # really should use QDesktopServices for this, but it appears to be broken.
        # url = QtCore.QUrl(self.experiment.lims_drawing_tool_url)
        # QtGui.QDesktopServices.openUrl(url)
        cmd = config.browser_command.format(url=self.experiment.lims_drawing_tool_url)
        subprocess.Popen(cmd, shell=True)

    def edit_pipettes_yml(self):
        pip_file = self.experiment.pipette_file
        if pip_file is None:
            raise Exception("No pipettes.yml file for this experiment.")
        cmd = config.editor_command.format(file=pip_file)
        subprocess.Popen(cmd, shell=True)


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
