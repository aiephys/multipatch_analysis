from sqlalchemy.orm import relationship
from . import make_table
from .dataset import PulseResponse, Baseline


__all__ = ['PulseResponseFit', 'PulseResponseStrength', 'BaselineResponseStrength']


PulseResponseFit = make_table(
    name='pulse_response_fit',
    comment="Curve fits to individual synaptic responses.",
    columns=[
        ('pulse_response_id', 'pulse_response.id', '', {'index': True, 'unique': True}),
        
        ('fit_amp', 'float', 'Fit amplitude of the response to this stimulus', {'index': True}),
        ('fit_latency', 'float', 'Fit latency of the response to this stimulus', {'index': True}),
        ('fit_yoffset', 'float', 'Fit y offset of the response to this stimulus'),
        ('fit_rise_time', 'float', 'Fit rise time of the response to this stimulus'),
        ('fit_decay_tau', 'float', 'Fit decay tau of the response to this stimulus'),
        ('fit_exp_amp', 'float', 'Fit exponential amplitude of the baseline before this stimulus'),
        ('fit_nrmse', 'float', 'Normalized RMS error of the fit'),

        ('baseline_fit_amp', 'float', 'Fit amplitude of the baseline before the stimulus', {'index': True}),
        ('baseline_fit_latency', 'float', 'Fit latency of the baseline before the stimulus', {'index': True}),
        ('baseline_fit_yoffset', 'float', 'Fit y offset of the baseline before the stimulus'),
        ('baseline_fit_rise_time', 'float', 'Fit rise time of the baseline before the stimulus'),
        ('baseline_fit_decay_tau', 'float', 'Fit decay tau of the baseline before the stimulus'),
        ('baseline_fit_exp_amp', 'float', 'Fit exponential amplitude of the baseline before this stimulus'),
        ('baseline_fit_nrmse', 'float', 'Normalized RMS error of the fit'),
    ]
)

PulseResponseStrength = make_table(
    name='pulse_response_strength',
    comment="Measurements of membrane potential or current deflection following each evoked presynaptic spike.",
    columns=[
        ('pulse_response_id', 'pulse_response.id', '', {'index': True, 'unique': True}),
        
        ('pos_amp', 'float', 'max-median offset from baseline to pulse response window'),
        ('neg_amp', 'float', 'min-median offset from baseline to pulse response window'),
        ('pos_dec_amp', 'float', 'max-median offset from baseline to pulse response window from devonvolved trace'),
        ('neg_dec_amp', 'float', 'min-median offset from baseline to pulse response window from deconvolved trace'),
        ('pos_dec_latency', 'float', 'duration (seconds) from presynaptic spike max dv/dt until the sample measured in pos_dec_amp'),
        ('neg_dec_latency', 'float', 'duration (seconds) from presynaptic spike max dv/dt until the sample measured in neg_dec_amp'),
        ('crosstalk', 'float', 'trace difference immediately before and after onset of presynaptic stimulus pulse'),
    ]
)

BaselineResponseStrength = make_table(
    name='baseline_response_strength',
    comment="Measurements of membrane potential or current deflection in the absence of presynaptic spikes (provides a measurement of background noise to compare to pulse_response_strength).",
    columns=[
        ('baseline_id', 'baseline.id', '', {'index': True, 'unique': True}),
        ('pos_amp', 'float', 'max-median offset from baseline to pulse response window'),
        ('neg_amp', 'float', 'min-median offset from baseline to pulse response window'),
        ('pos_dec_amp', 'float', 'max-median offset from baseline to pulse response window from devonvolved trace'),
        ('neg_dec_amp', 'float', 'min-median offset from baseline to pulse response window from deconvolved trace'),
        ('pos_dec_latency', 'float', 'duration (seconds) from presynaptic spike max dv/dt until the sample measured in pos_dec_amp'),
        ('neg_dec_latency', 'float', 'duration (seconds) from presynaptic spike max dv/dt until the sample measured in neg_dec_amp'),
        ('crosstalk', 'float', 'trace difference immediately before and after onset of presynaptic stimulus pulse'),
    ]
)

PulseResponse.pulse_response_fit = relationship(PulseResponseFit, back_populates="pulse_response", cascade="delete", single_parent=True, uselist=False)
PulseResponseFit.pulse_response = relationship(PulseResponse, back_populates="pulse_response_fit", single_parent=True)

PulseResponse.pulse_response_strength = relationship(PulseResponseStrength, back_populates="pulse_response", cascade="delete", single_parent=True, uselist=False)
PulseResponseStrength.pulse_response = relationship(PulseResponse, back_populates="pulse_response_strength", single_parent=True)

Baseline.baseline_response_strength = relationship(BaselineResponseStrength, back_populates="baseline", cascade="delete", single_parent=True, uselist=False)
BaselineResponseStrength.baseline = relationship(Baseline, back_populates="baseline_response_strength", single_parent=True)
