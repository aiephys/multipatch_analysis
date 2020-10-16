from sqlalchemy.orm import relationship, deferred, sessionmaker, aliased
from ... import config
from . import make_table
from .slice import Slice
from .experiment import Cell, Experiment

__all__ = ['CorticalCellLocation', 'CorticalSite']

CorticalCellLocation = make_table(
    name='cortical_cell_location',
    comment='Each row holds location information for a single cortical cell.',
    columns=[
        ('cell_id', 'cell.id', 'ID of the cell these locations apply to.', {'index':True}),
        ('cortical_site_id', 'cortical_site.id', 'ID of the site location measurements fit within.', {'index':True}),
        ('layer', 'str', 'Name of the layer the cell is in.'),
        ('distance_to_pia', 'float', 'The distance from the cell to the pial surface in m.'),
        ('distance_to_wm', 'float', 'The distance from the cell to the white matter in m.'),
        ('fractional_depth', 'float', 'The cortical depth of the cell where pia is 0 and wm is 1.'),
        ('layer_depth', 'float', 'Absolute depth within the layer in m.'),
        ('fractional_layer_depth', 'float', 'Fractional depth within the cells layer.'),
        ('position', 'object', '2D array, position of cell in slice image coordinates (in m)'),
        ])
        
Cell.cortical_location = relationship(CorticalCellLocation, back_populates="cell", cascade="delete", single_parent=True, uselist=False)
CorticalCellLocation.cell = relationship(Cell, back_populates="cortical_location", single_parent=True)

CorticalSite = make_table(
    name='cortical_site',
    comment='Each row holds measurements about one cortical site in a slice.',
    columns=[
        ('slice_id', 'slice.id', 'ID of the slice this site belongs to.', {'index':True}),
        ('experiment_id', 'experiment.id', 'ID to the experiment this site is part of', {'index':True}),
        ('pia_to_wm_distance', 'float', 'The distance (in m) from the pia to the white matter.'),
        ('pia_position', 'object', '3D location where the pia was marked in the arbitrary coordinate system of the experiment'),
        ('wm_position', 'object', '3D location where the wm was marked in the arbitrary coordinate system of the experiment'),
        ('layer_boundaries', 'object', 'Dictionary with fractional layer boundaries appropriate for the site.'),
        ('brain_region', 'str', 'The name of the brain region for the site.')
        ])

Slice.cortical_sites = relationship(CorticalSite, back_populates="slice", cascade="save-update,merge,delete")
CorticalSite.slice = relationship(Slice, back_populates='cortical_sites')

CorticalCellLocation.cortical_site = relationship(CorticalSite, back_populates='cell_locations')
CorticalSite.cell_locations = relationship(CorticalCellLocation, back_populates='cortical_site', cascade='delete')

CorticalSite.experiment = relationship(Experiment, back_populates='cortical_sites', single_parent=True)
Experiment.cortical_sites = relationship(CorticalSite, back_populates='experiment', cascade="save-update,merge,delete")