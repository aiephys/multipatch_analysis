from . import make_table
from sqlalchemy.orm import relationship
from .dataset import Recording
from .experiment import Cell

OpticalStimPulse = make_table(
    name='optical_stim_pulse',
    comment='A focal optical pulse stimulus',
    columns=[
        ('recording_id', 'recording.id', '', {'index': True}), ## this makes the assumption that the light is being recorded by a photodiode/pmt/other sensor
        ('pulse_number', 'int', 'The ordinal position of this pulse within a train of pulses.', {'index': True}),
        ('onset_time', 'float', 'The starting time of the pulse, relative to the beginning of the recording'),
        ('amplitude', 'float', 'Amplitude of the pulse'),
        ('duration', 'float', 'Length of the pulse in seconds'),
        ('wavelength', 'float', 'Wavelength of the light stimulus'),
        ('light_source', 'str', 'Name of the source of the light stimulus (LED, CoherentLaser, etc)'),
        ('cell_id', 'cell.id', 'Cell that was targeted by this stimulus, if any.'),
        ('position', 'object', '3D location of this stimulation in the arbitrary coordinate system of the experiment')
        ])

Recording.optical_stim_pulses = relationship(OpticalStimPulse, order_by=OpticalStimPulse.id, back_populates="recording", cascade='save-update,merge,delete')
OpticalStimPulse.recording = relationship(Recording, back_populates="optical_stim_pulses")
Cell.optical_stim_pulses = relationship(OpticalStimPulse, order_by=OpticalStimPulse.id, back_populates='cell', cascade='save-update,merge,delete')
OpticalStimPulse.cell = relationship(Cell, back_populates='optical_stim_pulses')