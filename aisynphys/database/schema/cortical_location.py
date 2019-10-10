from sqlalchemy.orm import relationship, deferred, sessionmaker, aliased
from ... import config
from . import make_table
from .slice import Slice
from .experiment import Cell, Experiment

__all__ = ['CellLocation', 'CorticalSite']

CellLocation = make_table(
    name='cell_location',
    comment='Each row holds location information for a single cell.',
    columns=[
        ('cell_id', 'cell.id', 'ID of the cell these locations apply to.', {'index':True}),
        ('cortical_site_id', 'cortical_site.id', 'ID of the site location measurements fit within.', {'index':True}),
        ('layer', 'str', 'Name of the layer the cell is in.'),
        ('distance_to_pia', 'float', 'The distance from the cell to the pial surface in m.'),
        ('distance_to_wm', 'float', 'The distance from the cell to the white matter in m.'),
        ('fractional_depth', 'float', 'The cortical depth of the cell where pia is 0 and wm is 1.')
        ])
        
Cell.location = relationship(CellLocation, back_populates="cell", cascade="delete", single_parent=True, uselist=False)
CellLocation.cell = relationship(Cell, back_populates="location", single_parent=True)

CorticalSite = make_table(
    name='cortical_site',
    comment='Each row holds measurements about one cortical site in a slice.',
    columns=[
        ('slice_id', 'slice.id', 'ID of the slice this site belongs to.', {'index':True}),
        ('experiment_id', 'experiment.id', 'ID to the experiment this site is part of', {'index':True}),
        ('pia_to_wm_distance', 'float', 'The distance (in m) from the pia to the white matter.'),
        ('pia_position', 'object', '3D location where the pia was marked in the arbitrary coordinate system of the experiment'),
        ('wm_position', 'object', '3D location where the wm was marked in the arbitrary coordinate system of the experiment'),
        ('L1_L23_boundary', 'float', 'Distance from the pia to the boundary between L1 and L2/3 as a fraction of cortex.'),
        ('L23_L4_boundary', 'float', 'Distance from the pia to the boundary between L2/3 and L4 as a fraction of cortex.'),
        ('L4_L5_boundary', 'float', 'Distance from the pia to the boundary between L4 and L5 as a fraction of cortex.'),
        ('L5_L6_boundary', 'float', 'Distance from the pia to the boundary between L5 and L6 as a fraction of cortex.'),
        ])

Slice.site = relationship(CorticalSite, back_populates="slice", cascade="delete", single_parent=True)
CorticalSite.slice = relationship(Slice, back_populates='site', single_parent=True)
CellLocation.site = relationship(CorticalSite, back_populates='cell')
CorticalSite.cell = relationship(CellLocation, back_populates='site', cascade='delete')
CorticalSite.experiment = relationship(Experiment, back_populates='site', single_parent=True)
Experiment.site = relationship(CorticalSite, back_populates='experiment', cascade="delete", single_parent=True)