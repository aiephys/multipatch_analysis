import numpy as np


class Pair(object):
    """Represents two cells that were recorded simultaneously; one putatively presynaptic
    and the other postsynaptic.
    
    Note: this is an _ordered_ pair; Pair(A, B) is not the same as Pair(B, A).
    
    Parameters
    ----------
    experiment : Experiment instance
        The experiment to which this pair belongs
    pre_cell : Cell instance
        Presynaptic cell
    post_cell : Cell instance
        Postsynaptic cell
    has_synapse : bool | None
        Whether a chemical synapse connects pre_cell to post_cell
    has_electrical : bool | None
        Whether an electrical synapse connects pre_cell to post_cell
    """
    def __init__(self, experiment, pre_cell, post_cell, has_synapse=None, synapse_type=None, has_electrical=None):
        self.experiment = experiment
        self.pre_cell = pre_cell
        self.post_cell = post_cell
        self.has_synapse = has_synapse
        self.has_electrical = has_electrical

    @property
    def distance(self):
        """Disance between cell positions
        """
        p1, p2 = self.pre_cell.position, self.post_cell.position
        if None in [p1, p2]:
            return None
        else:
            return np.linalg.norm(np.array(p1) - np.array(p2))
