from neuroanalysis.ui.nwb_viewer import MiesNwbViewer
from neuroanalysis.miesnwb import MiesNwb

from .multipatch_view import MultipatchMatrixView
from aisynphys.data import MultiPatchDataset
from .pair_view import PairView


class MultipatchNwbViewer(MiesNwbViewer):
    def create_views(self):
        MiesNwbViewer.create_views(self)
        
        add_views = [
            ('Matrix', MultipatchMatrixView(self)),
            ('Pair', PairView(self)),
        ]
        for name, view in add_views:
            self.tabs.addTab(view, name)
    
    def load_nwb(self, filename):
        nwb = MultiPatchDataset(filename)
        self.set_nwb(nwb)
        return nwb
