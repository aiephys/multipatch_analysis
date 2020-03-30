from sqlalchemy.orm import relationship
from . import make_table
from .experiment import Pair


__all__ = ['GapJunction']


GapJunction = make_table(
    name='gap_junction',
    comment= "Describes the presence of gap junction.",
    columns=[
        ('pair_id', 'pair.id', 'The ID of the entry in the pair table to which these results apply', {'index': True}),
        ('corr_coeff_pulse', 'float', 'The Pearson correlation coefficient of pre- and post-synaptic long pulse', {'index': True}),
        ('corr_coeff_noise', 'float', 'The Pearson correlation coefficient of pre- and post-synaptic background', {'index': True}),
        ('p_val_pulse', 'float', 'The Pearson p-value of pre- and post-synaptic long pulse', {'index': True}),
        ('p_val_noise', 'float', 'The Pearson p-value of pre- and post-synaptic background', {'index': True}),
        ('coupling_coeff_pulse', 'float', 'The coupling coefficient of pre- and post-synaptic long pulse', {'index': True}),
        ('coupling_coeff_noise', 'float', 'The coupling coefficient of pre- and post-synaptic background', {'index': True}),
    ]
)

Pair.gap_junction = relationship(GapJunction, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
GapJunction.pair = relationship(Pair, back_populates="gap_junction", single_parent=True)
