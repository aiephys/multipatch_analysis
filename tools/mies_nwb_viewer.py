import os, sys
import acq4
from acq4.pyqtgraph.Qt import QtCore, QtGui
import acq4.pyqtgraph.console

from multipatch_analysis.ui.multipatch_nwb_viewer import MultipatchNwbViewer


if __name__ == '__main__':
    # open a console for debugging
    message = """Console variables:

    analyzer : The currently active analyzer
    sweeps : A list of the currently selected sweeps
    view : The MiesNwbViewer instance
    nwb : The currently loaded NWB file
    man : ACQ4 manager
    man.currentFile : currently selected file in ACQ4 DataManager module

    """
    console = acq4.pyqtgraph.console.ConsoleWidget(text=message)
    console.catchNextException()
    console.show()

    # start up ACQ4 data manager
    m = acq4.Manager.Manager(argv=['-D', '-n', '-m', 'Data Manager'])
    dm = m.getModule('Data Manager')

    # create NWB viewer
    v = MultipatchNwbViewer()
    v.show()
    v.setWindowTitle('Multipatch NWB Viewer')

    # set up a convenient function for loading nwb from filename
    nwb = None
    def load_nwb(filename):
        global nwb
        print("Loading %s" % filename)
        nwb = v.load_nwb(filename)
        console.localNamespace['nwb'] = nwb

    # add a button to load nwb from file selected in data manager
    def load_from_dm():
        with acq4.pyqtgraph.BusyCursor():
            load_nwb(m.currentFile.name())

    btn = acq4.pyqtgraph.Qt.QtGui.QPushButton('load from data manager')
    v.vsplit.insertWidget(0, btn)
    btn.clicked.connect(load_from_dm)

    # make a few variables available from the console
    console.localNamespace.update({'man': m, 'view': v, 'nwb': None, 'sweeps': [], 'analyzer': None})

    # make selected sweeps / analyzer available as well
    def update_namespace():
        console.localNamespace.update({'sweeps': v.selected_sweeps(), 'analyzer': v.selected_analyzer()})

    v.analyzer_changed.connect(update_namespace)
    v.explorer.selection_changed.connect(update_namespace)

    # Set up code reloading shortcut
    def reload_all():
        acq4.pyqtgraph.reload.reloadAll(verbose=True)

    reload_shortcut = QtGui.QShortcut(QtGui.QKeySequence('Ctrl+r'), v)
    reload_shortcut.setContext(QtCore.Qt.ApplicationShortcut)
    reload_shortcut.activated.connect(reload_all)

    # load file or set base directory from argv
    for arg in sys.argv[1:]:
        if not os.path.exists(arg):
            from multipatch_analysis.experiment_list import cached_experiments
            expts = cached_experiments()
            try:
                expt = expts[arg]
            except Exception:
                print("Could not find experiment %s" % arg)
                sys.exit(-1)
            arg = expt.nwb_file

        if os.path.isdir(arg):
            m.setBaseDir(arg)
        elif os.path.isfile(arg):
            load_nwb(arg)

    # start Qt event loop if this is not an interactive python session
    if sys.flags.interactive == 0:
        acq4.pyqtgraph.QtGui.QApplication.exec_()
