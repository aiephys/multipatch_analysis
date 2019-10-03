from sqlalchemy.orm import relationship
from . import make_table
from .experiment import Pair


__all__ = ['Synapse', 'AvgResponseFit']


Synapse = make_table(
    name='synapse',
    comment="Chemical synapse properties",
    columns=[
        ('pair_id', 'pair.id', 'The ID of the entry in the pair table to which these results apply', {'index': True}),
        ('synapse_type', 'str', '"ex" or "in" indicating whether the synapse is excitatory or inhibitory', {'index': True}),
        ('latency', 'float', 'Latency in seconds from spike max slope until synaptic response onset.', {'index': True}),
        ('psp_amplitude', 'float', 'Amplitude of resting-state PSPs in Volts.'),
        ('psp_rise_time', 'float', 'Rise time in seconds measured from averaged PSPs.'),
        ('psp_decay_tau', 'float', 'decay time constant in seconds measured from averaged PSPs.'),
        ('psc_amplitude', 'float', 'Amplitude of resting-state PSCs in Amperes.'),
        ('psc_rise_time', 'float', 'Rise time in seconds measured from averaged PSCs.'),
        ('psc_decay_tau', 'float', 'decay time constant in seconds measured from averaged PSCs.'),
    ]
)


AvgResponseFit = make_table(
    name='avg_response_fit',
    comment="Fit to average post synaptic response for a given pair. Each pair may have fits for VC and IC recordings, held at -70 and -55 mV.",
    columns=[
        ('pair_id', 'pair.id', 'The ID of the entry in the pair table to which these results apply', {'index': True}),
        ('clamp_mode', 'str', 'The clamp mode "ic" or "vc"', {'index': True}),
        ('holding', 'float', 'The holding potential -70 or -55', {'index': True}),
        ('fit_xoffset', 'float', 'Fit time from max slope of the presynaptic spike until onset of the synaptic response (seconds)'),
        ('fit_yoffset', 'float', 'Fit constant y-offset (amps or volts)'),
        ('fit_amp', 'float', 'Fit synaptic response amplitude (amps or volts)'),
        ('fit_rise_time', 'float', 'Fit rise time (seconds) from response onset until peak'),
        ('fit_rise_power', 'float', 'Fit rise exponent (usually fixed at 2)'),
        ('fit_decay_tau', 'float', 'Fit exponential decay time constant (seconds)'),
        ('fit_exp_amp', 'float', 'Fit baseline exponental amplitude (amps or volts)'),
        ('fit_exp_tau', 'float', 'Fit baseline exponental decay time constant (seconds)'),
        ('nrmse', 'float', 'Normalized RMS error of the fit residual'),
        ('initial_xoffset', 'float', 'Initial latency supplied to fitting algorithm'),
        ('manual_qc_pass', 'bool', 'If true, this fit passes manual verification QC'),
        ('avg_data', 'array', 'Averaged PSP/PSC that was fit.', {'deferred': True}),
        ('avg_data_start_time', 'float', 'Starting time of avg_data, relative to the presynaptic spike'),
        ('n_averaged_responses', 'int', 'Number of postsynaptic responses that were averaged in avg_data'),
        ('avg_baseline_noise', 'float', 'Standard deviation of avg_data before the presynaptic stimulus'),
    ]
)


Pair.synapse = relationship(Synapse, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
Synapse.pair = relationship(Pair, back_populates="synapse", single_parent=True)

Pair.avg_response_fits = relationship(AvgResponseFit, back_populates="pair", cascade="delete", single_parent=True, uselist=True)
AvgResponseFit.pair = relationship(Pair, back_populates="avg_response_fits", single_parent=True)


