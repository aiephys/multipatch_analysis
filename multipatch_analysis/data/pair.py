

class Pair(object):
    """Represents two cells that were recorded simultaneously; one putatively presynaptic
    and the other postsynaptic.
    
    Note: this is an _ordered_ pair; Pair(A, B) is not the same as Pair(B, A).
    
    Parameters
    ----------
    pre_cell : Cell instance
        Presynaptic cell
    post_cell : Cell instance
        Postsynaptic cell
    synapse : bool | None
        Whether a chemical synapse connects pre_cell to post_cell
    synapse_type : str | None
        'ex' or 'in' indicating whether the synapse is excitatory or inhibitory
    electrical : bool | None
        Whether an electrical synapse connects pre_cell to post_cell
    """
    def __init__(self, pre_cell, post_cell, synapse=None, synapse_sign=None, electrical=None):
        self.pre_cell = pre_cell
        self.post_cell = post_cell
        self.synapse = synapse
        self.synapse_sign = synapse_sign
        self.electrical = electrical

    @property
    def distance(self):
        p1, p2 = self.pre_cell.position, self.post_cell.position
        if None in [p1, p2]:
            return None
        else:
            return np.linalg.norm(np.array(p1) - np.array(p2))
