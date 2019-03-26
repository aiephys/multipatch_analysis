from sqlalchemy.orm import relationship
from .database import TableGroup, make_table, Column, Integer, Boolean, ForeignKey
from .experiment import Cell


__all__ = ['morphology_tables', 'Morphology']

Morphology = make_table(name='morphology', comment="Describes morphological properties of cells.", columns=[
    ('cell_id', 'cell.id', 'The ID of the cell described by each record', {'index': True, 'unique': True}),
    ('pyramidal', 'bool', 'Whether the experimenter labeled this cell as pyramidal', {'index': True}),
])

Cell.morphology = relationship(Morphology, back_populates="cell", cascade="delete", single_parent=True, uselist=False)
Morphology.cell = relationship(Cell, back_populates="morphology", single_parent=True)

morphology_tables = TableGroup([Morphology])
