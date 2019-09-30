from sqlalchemy.orm import relationship
from . import make_table
from .experiment import Pair

__all__ = ['RestingStateFit']


RestingStateFit = make_table(
    name='resting_state_fit',
    comment="""Contains curve fits to averages of "resting state" synaptic responses.""",
    columns=[
        ('pair_id', 'pair.id', 'The ID of the entry in the pair table to which these results apply', {'index': True}),

        # current clamp
        ('ic_amp', 'float', 'fit amplitude of current clamp average first pulses'),
        ('ic_latency', 'float', 'fit time elapsed since the time of presynaptic spike (max dv/dt) of current clamp data'),
        ('ic_rise_time', 'float', 'fit rise time of psp of current clamp data'),
        ('ic_decay_tau', 'float', 'fit decay of psp of current clamp data'),
        ('ic_exp_amp', 'float', 'fit amplitude of exponential baseline'),
        ('ic_exp_tau', 'float', 'fit tau of exponential baseline'),
        ('ic_avg_data', 'array', 'array of the data voltage waveform used in fitting'),
        ('ic_avg_data_start_time', 'float', 'time value of the first sample in ic_avg_data, relative to the presynaptic spike'),
        ('ic_pulse_ids', 'array', 'data base pulse ids included in the current clamp fit'),
        ('ic_nrmse', 'float', 'error of fit of current clamp fit'),

        # voltage clamp
        ('vc_amp', 'float', 'fit amplitude of voltage clamp average first pulses'),
        ('vc_latency', 'float', 'fit time elapsed since the time of presynaptic spike (max dv/dt) of voltage clamp data'),
        ('vc_rise_time', 'float', 'fit rise time of psp measured in voltage clamp'),
        ('vc_decay_tau', 'float', 'fit decay of psp measured in voltage clamp'),
        ('vc_exp_amp', 'float', 'fit amplitude of exponential baseline'),
        ('vc_exp_tau', 'float', 'fit tau of exponential baseline'),
        ('vc_avg_data', 'array', 'array of the data current waveform used in fitting'),
        ('vc_avg_data_start_time', 'float', 'time value of the first sample in vc_avg_data, relative to the presynaptic spike'),
        ('vc_pulse_ids', 'array', 'data base pulse ids included in the voltage clamp fit'),
        ('vc_nrmse', 'float', 'error of fit of voltage clamp fit'),
    ]
)

Pair.resting_state_fit = relationship(RestingStateFit, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
RestingStateFit.pair = relationship(Pair, back_populates="resting_state_fit", single_parent=True)
