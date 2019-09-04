from sqlalchemy.orm import relationship
from . import make_table
from .pulse_response import PulseResponseStrength, BaselineResponseStrength
from .experiment import Pair


__all__ = ['SynapsePrediction']


SynapsePrediction = make_table(
    name='synapse_prediction',
    comment= "Unbiased metrics used for automated synapse detection.",
    columns=[
        ('pair_id', 'pair.id', 'The ID of the entry in the pair table to which these results apply', {'index': True}),
        ('synapse_type', 'str', 'String "ex" or "in", indicating whether this analysis chose to treat the pair as excitatory or inhibitory'),

        # current clamp metrics
        ('ic_n_samples', 'int', "Number of samples (pulse responses) that were pooled from current clamp recordings"),
        ('ic_crosstalk_mean', 'float'),
        ('ic_base_crosstalk_mean', 'float'),
        # amplitude,
        ('ic_amp_mean', 'float'),
        ('ic_amp_stdev', 'float'),
        ('ic_base_amp_mean', 'float'),
        ('ic_base_amp_stdev', 'float'),
        ('ic_amp_ttest', 'float'),
        ('ic_amp_ks2samp', 'float'),
        # deconvolved amplitide
        ('ic_deconv_amp_mean', 'float'),
        ('ic_deconv_amp_stdev', 'float'),
        ('ic_base_deconv_amp_mean', 'float'),
        ('ic_base_deconv_amp_stdev', 'float'),
        ('ic_deconv_amp_ttest', 'float'),
        ('ic_deconv_amp_ks2samp', 'float'),
        # latency
        ('ic_latency_mean', 'float'),
        ('ic_latency_stdev', 'float'),
        ('ic_base_latency_mean', 'float'),
        ('ic_base_latency_stdev', 'float'),
        ('ic_latency_ttest', 'float'),
        ('ic_latency_ks2samp', 'float'),
        
        # voltage clamp metrics
        ('vc_n_samples', 'int'),
        ('vc_crosstalk_mean', 'float'),
        ('vc_base_crosstalk_mean', 'float'),
        # amplitude,
        ('vc_amp_mean', 'float'),
        ('vc_amp_stdev', 'float'),
        ('vc_base_amp_mean', 'float'),
        ('vc_base_amp_stdev', 'float'),
        ('vc_amp_ttest', 'float'),
        ('vc_amp_ks2samp', 'float'),
        # deconvolved amplitide
        ('vc_deconv_amp_mean', 'float'),
        ('vc_deconv_amp_stdev', 'float'),
        ('vc_base_deconv_amp_mean', 'float'),
        ('vc_base_deconv_amp_stdev', 'float'),
        ('vc_deconv_amp_ttest', 'float'),
        ('vc_deconv_amp_ks2samp', 'float'),
        # latency
        ('vc_latency_mean', 'float'),
        ('vc_latency_stdev', 'float'),
        ('vc_base_latency_mean', 'float'),
        ('vc_base_latency_stdev', 'float'),
        ('vc_latency_ttest', 'float'),
        ('vc_latency_ks2samp', 'float'),
        
        # Average pulse responses
        ('ic_average_response', 'array'),
        ('ic_average_response_t0', 'float'),
        ('ic_average_base_stdev', 'float'),
        ('vc_average_response', 'array'),
        ('vc_average_response_t0', 'float'),
        ('vc_average_base_stdev', 'float'),

        # PSP fit parameters
        ('ic_fit_amp', 'float'),
        ('ic_fit_xoffset', 'float'),
        ('ic_fit_yoffset', 'float'),
        ('ic_fit_rise_time', 'float'),
        ('ic_fit_rise_power', 'float'),
        ('ic_fit_decay_tau', 'float'),
        ('ic_fit_exp_amp', 'float'),
        ('ic_fit_exp_tau', 'float'),
        ('ic_fit_nrmse', 'float'),

        ('vc_fit_amp', 'float'),
        ('vc_fit_xoffset', 'float'),
        ('vc_fit_yoffset', 'float'),
        ('vc_fit_rise_time', 'float'),
        ('vc_fit_rise_power', 'float'),
        ('vc_fit_decay_tau', 'float'),
        ('vc_fit_exp_amp', 'float'),
        ('vc_fit_exp_tau', 'float'),
        ('vc_fit_nrmse', 'float'),
    ]
)

Pair.synapse_prediction = relationship(SynapsePrediction, back_populates="pair", cascade="delete", single_parent=True, uselist=False)
SynapsePrediction.pair = relationship(Pair, back_populates="synapse_prediction", single_parent=True)

