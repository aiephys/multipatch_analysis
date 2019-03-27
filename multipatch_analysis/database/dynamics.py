from sqlalchemy.orm import relationship
from .database import make_table, TableGroup
from .experiment import Pair


__all__ = ['dynamics_tables', 'Dynamics']


Dynamics = make_table(
    name='dynamics',
    comment="Describes short term dynamics of synaptic connections.",
    columns=[
        ('pair_id', 'pair.id', 'The ID of the cell pair described by each record', {'index': True, 'unique': True}),
        ('pulse_ratio_8_1_50Hz', 'float', '8:1 pulse ratio for 50Hz induction', {'index': True}),
        ('pulse_ratio_2_1_50Hz', 'float', '2:1 pulse ratio for 50Hz induction', {'index': True}),
    ]
)

Pair.dynamics = relationship(Dynamics, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
Dynamics.pair = relationship(Pair, back_populates="dynamics", single_parent=True)


dynamics_tables = TableGroup([Dynamics])
