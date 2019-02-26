from sqlalchemy.orm import relationship
from .database import TableGroup
from .experiment import Cell


class MorphologyTableGroup(TableGroup):
    schemas = {
        'morphology': [
            """Describes morphological properties of cells.
            """,
            ('cell_id', 'cell.id', 'The ID of the cell described by each record', {'index': True, 'unique': True}),
            ('pyramidal', 'bool', 'Whether the experimenter labeled this cell as pyramidal', {'index': True}),
            # TODO: import more features from LIMS
        ],
    }

    def create_mappings(self):
        TableGroup.create_mappings(self)
        
        Morphology = self['morphology']
        
        Cell.morphology = relationship(Morphology, back_populates="cell", cascade="delete", single_parent=True, uselist=False)
        Morphology.cell = relationship(Cell, back_populates="morphology", single_parent=True)


morphology_tables = MorphologyTableGroup()


def init_tables():
    global Morphology
    morphology_tables.create_tables()
    Morphology = morphology_tables['morphology']


# create tables in database and add global variables for ORM classes
init_tables()
