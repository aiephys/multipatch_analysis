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
