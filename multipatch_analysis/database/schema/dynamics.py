from sqlalchemy.orm import relationship
from . import make_table
from .experiment import Pair


__all__ = ['Dynamics']


Dynamics = make_table(
    name='dynamics',
    comment="Describes short term dynamics of synaptic connections.",
    columns=[
        ('pair_id', 'pair.id', 'The ID of the cell pair described by each record', {'index': True, 'unique': True}),
        ('paired_pulse_ratio_50hz', 'float', 'The mean ratio of 2nd / 1st pulse amplitudes for 50Hz pulse trains.', {'index': True}),
        ('stp_initial_50hz', 'float', 'The mean relative change from 1st to 2nd pulse for 50Hz pulse trains', {'index': True}),
        ('stp_initial_50hz_n', 'float', 'Number of samples represented in stp_initial_50Hz', {'index': True}),
        ('stp_initial_50hz_std', 'float', 'Standard deviation of samples represented in stp_initial_50Hz', {'index': True}),
        ('stp_induction_50hz', 'float', 'The mean relative change from 1st to 5th-8th pulses for 50Hz pulse trains', {'index': True}),
        ('stp_induction_50hz_n', 'float', 'Number of samples represented in stp_induction_50Hz', {'index': True}),
        ('stp_induction_50hz_std', 'float', 'Standard deviation of samples represented in stp_induction_50Hz', {'index': True}),
        ('stp_recovery_250ms', 'float', 'The mean relative change from 1st-4th to 9th-12th pulses for pulse trains with a 250 ms recovery period', {'index': True}),
        ('stp_recovery_250ms_n', 'float', 'Number of samples represented in stp_recovery_250ms', {'index': True}),
        ('stp_recovery_250ms_std', 'float', 'Standard deviation of samples represented in stp_recovery_250ms', {'index': True}),
        ('pulse_amp_90th_percentile', 'float', 'The 90th-percentile largest pulse amplitude, used to normalize change values in this table', {}),
        ('noise_amp_90th_percentile', 'float', 'The 90th-percentile largest amplitude measured from background noise, used for comparison to pulse_amp_90th_percentile', {}),
    ]
)

Pair.dynamics = relationship(Dynamics, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
Dynamics.pair = relationship(Pair, back_populates="dynamics", single_parent=True)
