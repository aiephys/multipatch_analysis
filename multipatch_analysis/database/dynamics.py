from sqlalchemy.orm import relationship
from .database import TableGroup
from .experiment import Pair


__all__ = ['dynamics_tables', 'Dynamics']


class DynamicsTableGroup(TableGroup):
    schemas = {
        'dynamics': [
            {'comment': """Describes short term dynamics of connections.
            """},
            ('pair_id', 'pair.id', 'The ID of the cell pair described by each record', {'index': True, 'unique': True}),
            ('pulse_ratio_8_1_50Hz', 'float', '8:1 pulse ratio for 50Hz induction', {'index': True}),
            ('pulse_ratio_2_1_50Hz', 'float', '2:1 pulse ratio for 50Hz induction', {'index': True}),
        ],
    }

    def create_mappings(self):
        TableGroup.create_mappings(self)
        
        Dynamics = self['dynamics']
        
        Pair.dynamics = relationship(Dynamics, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
        Dynamics.pair = relationship(Pair, back_populates="dynamics", single_parent=True)


dynamics_tables = DynamicsTableGroup()
Dynamics = dynamics_tables['dynamics']
