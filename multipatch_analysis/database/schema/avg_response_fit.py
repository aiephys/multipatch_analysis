from sqlalchemy.orm import relationship
from . import make_table
from .dataset import AvgResponseFit

__all__ = ['AvgResponseFit']

AvgResponseFit = make_table(
    name='avg_response_fit',
    comment="Fit to average post synaptic response for a given pair",
    columns=[
        ('pair_id', 'pair.id', 'The ID of the entry in the pair table to which these results apply', {'index': True}),
        ('latency', 'float', 'Latency supplied to fitting algorithm'),
        ('fit_parameters', 'object', 'Dictionary of fit parameters'),
    ]

Pair.avg_response_fit = relationship(AvgResponseFit, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
AvgResponseFit.pair = relationship(Pair, back_populates="avg_response_fit", single_parent=True)