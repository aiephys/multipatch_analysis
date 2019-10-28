from sqlalchemy.orm import relationship
from . import make_table
from .experiment import Cell


__all__ = ['Morphology']

Morphology = make_table(name='morphology', comment="Describes morphological properties of cells.", columns=[
    ('cell_id', 'cell.id', 'The ID of the cell described by each record', {'index': True, 'unique': True}),
    ('pyramidal', 'bool', 'Indicates whether the experimenter labeled the cell as pyramidal. '
        'This call is based on the presence of a prominent apical dendrite seen in the fluorescent dye fill during experiment. '
        'The `dendrite_type` column is recommended as a more reliable indicator of excitatory morphology.', {'index': True}),
    ('cortical_layer', 'str', 'Cortical layer of cell defined by layer drawing annotation based on DAPI staining', {'index': True}),
    ('qual_morpho_type', 'str', 'Qualitative desctription of cell morphology', {'index': True}),
    ('dendrite_type', 'str', 'Dendrite type of cell (spiny, aspiny, sparsely spiny) determined from biocytin staining. '
        'Generally spiny cells are taken to be excitatory; aspiny or sparsely spiny cells are inhibitory.', {'index': True}),
    ('apical_trunc_distance', 'float', 'Measured distance to truncation of apical dendrite', {'index': True}),
    ('axon_trunc_distance', 'float', 'Measured distance to truncation of axon', {'index': True}),
    ('apical_truncation', 'str', 'Qualitative description of apical dendrite truncation', {'index': True}),
    ('axon_truncation', 'str', 'Qualitative description of axon truncation', {'index': True}),
    ('axon_origin', 'str', 'Origination of axon; soma, dendrite, etc', {'index': True}),
    ('morpho_db_hash', 'bigint', 'hash of Morphology team database for update tracking'),
])

Cell.morphology = relationship(Morphology, back_populates="cell", cascade="delete", single_parent=True, uselist=False)
Morphology.cell = relationship(Cell, back_populates="morphology", single_parent=True)
