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
        ('input_resistance', 'float', 'Input resistance of the cell (from response peaks, capturing properties of channels at baseline)', {'index': True}),
        ('input_resistance_ss', 'float', '(True) steady-state input resistance of the cell', {'index': True}),
        ('sag', 'float', 'Hyperpolarizing sag ratio (peak/steady-state), measured from ~ -100mV current injection', {'index': True}),
        ('tau', 'float', 'Membrane time constant', {'index': True}),
        ('sag_peak_t', 'float', 'Time of peak hyperpolarizing sag.', {'index': True}),
        ('sag_depol', 'float', 'Depolarizing sag ratio (peak/steady-state), measured from the largest subthreshold depolarizing input', {'index': True}),
        ('sag_peak_t_depol', 'float', 'Time of peak depolarizing sag.', {'index': True}),
        ('ap_upstroke_downstroke_ratio', 'float', 'The upstroke-downstroke ratio of the first spike', {'index': True}),
        ('ap_width', 'float', 'Spike width', {'index': True}),
        ('ap_upstroke', 'float', 'Spike upstroke rate', {'index': True}),
        ('ap_downstroke', 'float', 'Spike downstroke rate', {'index': True}),
        ('ap_threshold_v', 'float', 'Spike threshold voltage', {'index': True}),
        ('ap_peak_deltav', 'float', 'Spike peak voltage relative to threshold', {'index': True}),
        ('ap_fast_trough_deltav', 'float', 'AHP / fast trough voltage relative to threshold', {'index': True}),

        ('firing_rate_rheo', 'float', 'Mean firing rate for rheobase sweep', {'index': True}),
        ('latency_rheo', 'float', 'First spike latency for rheobase sweep', {'index': True}),
        ('firing_rate_40pa', 'float', 'Mean firing rate for +40pA sweep (relative to rheobase)', {'index': True}),
        ('latency_40pa', 'float', 'First spike latency for +40pA sweep (relative to rheobase)', {'index': True}),
        
        ('adaptation_index', 'float', 'Adaptation index (ratio of consecutive ISIs), averaged across sweeps', {'index': True}),
        ('isi_cv', 'float', 'Coefficient of variation of ISI distribution, averaged across sweeps', {'index': True}),
        
        ('chirp_peak_freq', 'float', 'Frequency at which the chirp response peaks', {'index': True}),
        ('chirp_3db_freq', 'float', 'Frequency at which the chirp response amplitude is 3 dB below the peak.', {'index': True}),
        ('chirp_peak_ratio', 'float', 'Ratio of chirp resonance peak amplitude to low-frequency response amplitude', {'index': True}),
        ('chirp_peak_impedance', 'float', 'Impedance at chirp resonance peak.', {'index': True}),
        ('chirp_sync_freq', 'float', 'Frequency at which the chirp phase response equals zero.', {'index': True}),
        ('chirp_inductive_phase', 'float', 'Integrated of chirp phase response where phase > 0 (below sync freq).', {'index': True}),
        
        ('isi_adapt_ratio', 'float', 'Ratio of ISI on 5th to 1st spike', {'index': True}),
        ('upstroke_adapt_ratio', 'float', 'Ratio of upstroke on 5th to 1st spike', {'index': True}),
        ('downstroke_adapt_ratio', 'float', 'Ratio of downstroke on 5th to 1st spike', {'index': True}),
        ('width_adapt_ratio', 'float', 'Ratio of spike width on 5th to 1st spike', {'index': True}),
        ('threshold_v_adapt_ratio', 'float', 'Ratio of spike threshold on 5th to 1st spike', {'index': True}),
    ]
)

Cell.intrinsic = relationship(Intrinsic, back_populates="cell", cascade="delete", single_parent=True, uselist=False)
Intrinsic.cell = relationship(Cell, back_populates="intrinsic", single_parent=True)