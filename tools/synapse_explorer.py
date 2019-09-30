import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
from neuroanalysis.ui.plot_grid import PlotGrid
from aisynphys.ui.ui import SynapseTreeWidget, ExperimentInfoWidget
from aisynphys.synaptic_dynamics import DynamicsAnalyzer


class SynapseExplorer(QtGui.QWidget):
    def __init__(self, expts, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.expts = expts

        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        
        self.hsplit = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.hsplit, 0, 0)
        
        self.vsplit = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.hsplit.addWidget(self.vsplit)
        
        self.syn_tree = SynapseTreeWidget(self.expts)
        self.vsplit.addWidget(self.syn_tree)
        self.syn_tree.itemSelectionChanged.connect(self.selection_changed)
        
        self.expt_info = ExperimentInfoWidget()
        self.vsplit.addWidget(self.expt_info)
        
        self.train_plots = PlotGrid()
        self.hsplit.addWidget(self.train_plots)
        
        self.analyzers = {}
        
    def selection_changed(self):
        with pg.BusyCursor():
            sel = self.syn_tree.selectedItems()[0]
            expt = sel.expt
            
            self.expt_info.set_experiment(expt)
            
            pre_cell = sel.cells[0].cell_id
            post_cell = sel.cells[1].cell_id
            
            key = (expt, pre_cell, post_cell)
            if key not in self.analyzers:
                self.analyzers[key] = DynamicsAnalyzer(*key)
            analyzer = self.analyzers[key]
            
            if len(analyzer.pulse_responses) == 0:
                raise Exception("No suitable data found for cell %d -> cell %d in expt %s" % (pre_cell, post_cell, expt))
            
            # Plot all individual and averaged train responses for all sets of stimulus parameters
            self.train_plots.clear()
            analyzer.plot_train_responses(plot_grid=self.train_plots)
        


if __name__ == '__main__':
    from aisynphys.experiment_list import cached_experiments
    import sys
    app = pg.mkQApp()

    all_expts = cached_experiments()
    
    exp = SynapseExplorer(all_expts)
    exp.show()


    if sys.flags.interactive == 0:
        app.exec_()
