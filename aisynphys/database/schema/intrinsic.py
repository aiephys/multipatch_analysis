from sqlalchemy.orm import relationship
from . import make_table
from .experiment import Cell


__all__ = ['Intrinsic']


Intrinsic = make_table(
    name='intrinsic',
    comment= "Describes the intrinsic properties of cells using 1 sec long current steps and chirps.",
    columns=[
        ('cell_id', 'cell.id', 'The ID of the entry in the cell table to which these results apply', {'index': True}),
        ('rheobase', 'float', 'Current at rheobase', {'index': True}),
        ('fi_slope', 'float', 'Slope of the current-spiking relationship', {'index': True}),
        ('input_resistance', 'float', 'Input resistance of the cell', {'index': True}),
        ('sag', 'float', 'Hyperpolarizing sag measured from ~ -100mV current injection', {'index': True}),
        ('adaptation_index', 'float', 'Average adaptation index', {'index': True}),
        ('upstroke_downstroke_ratio', 'float', 'The upstroke-downstroke ratio of the first spike', {'index': True}),
        ('width', 'float', 'Spike width', {'index': True}),
        ('upstroke', 'float', 'Spike upstroke rate', {'index': True}),
        ('downstroke', 'float', 'Spike downstroke rate', {'index': True}),
        ('peak_v', 'float', 'Spike peak voltage', {'index': True}),
        ('threshold_v', 'float', 'Spike threshold voltage', {'index': True}),
        ('fast_trough_v', 'float', 'AHP / fast trough voltage', {'index': True}),
        ('chirp_peak_freq', 'float', 'Frequency at which the chirp response peaks', {'index': True}),
        ('chirp_3db_freq', 'float', 'Frequency at which the chirp response amplitude is 3 dB below the peak.', {'index': True}),
        ('chirp_peak_ratio', 'float', 'Ratio of chirp resonance peak amplitude to low-frequency response amplitude', {'index': True}),
        
    ]
)

Cell.intrinsic = relationship(Intrinsic, back_populates="cell", cascade="delete", single_parent=True, uselist=False)
Intrinsic.cell = relationship(Cell, back_populates="intrinsic", single_parent=True)