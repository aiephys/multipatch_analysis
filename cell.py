from __future__ import print_function, division
import numpy as np

from constants import ALL_CRE_TYPES, ALL_LABELS


class Cell(object):
    def __init__(self, expt, cell_id):
        self.expt = expt
        self.cell_id = cell_id
        self.access_qc = None
        self.holding_qc = None
        self.spiking_qc = None
        self.labels = {}
        self.position = None

    @property
    def pass_qc(self):
        """True if cell passes QC.
        """
        if self.access_qc is True and self.holding_qc is True:
            return True
        elif self.access_qc is False or self.holding_qc is False:
            return False

        # None means cell is not present in ephys data
        return None

    @property
    def cre_type(self):
        """Cre type string for this cell.
        
        If the cell is reporter-negative then cre_type is 'unk'.
        If the cell has ambiguous or missing data then cre_type is None.
        """
        default = 'unknown'
        ct = None
        for label,pos in self.labels.items():
            if label in ALL_LABELS:
                continue
            if pos == '+':
                if ct not in (None, default):
                    raise Exception("%s has multiple labels!" % self)
                ct = label
            elif pos == '-':
                if ct is not None:
                    continue
                ct = default
        return ct

    @property
    def label_type(self):
        """fluorescent type string for this cell.
        
        If the cell is reporter-negative then cre_type is 'unk'.
        If the cell has ambiguous or missing data then cre_type is None.
        """
        default = 'unknown'
        ct = None
        for label,pos in self.labels.items():
            if label in ALL_CRE_TYPES or label == 'biocytin':
                continue
            if pos == '+':
                if ct not in (None, default):
                    raise Exception("%s has multiple labels!" % self)
                ct = label
            elif pos == '-':
                if ct is not None:
                    continue
                ct = default
        return ct

    @property
    def depth(self):
        """Depth of cell from the cut surface of the slice.
        """
        sd = self.expt.surface_depth
        p = self.position
        if None in (sd, p):
            return None
        return sd - p[2]

    def distance(self, cell):
        """Return distance between cells, or nan if positions are not defined.
        """
        p1 = self.position
        p2 = cell.position
        if p1 is None or p2 is None:
            return np.nan
        return ((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)**0.5

    def __repr__(self):
        return "<Cell %s:%d>" % (self.expt.expt_id, self.cell_id)