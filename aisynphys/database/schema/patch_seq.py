from sqlalchemy.orm import relationship
from . import make_table
from .experiment import Cell


__all__ = ['PatchSeq']

PatchSeq = make_table(
    name='patch_seq', 
    comment="Describes transcriptomic data of a cell obtained from patch_seq experiments", 
    columns=[
    # These values are provided at the time of the experiment
    ('cell_id', 'cell.id', 'The ID of the cell described by each record', {'index': True, 'unique': True}),
    ('tube_id', 'str', 'Patched Cell Container ID used for RNA analysis', {'index': True}),
    ('nucleus', 'bool', 'Whether the nucleus was recovered from the cell', {'index': True}),
    ('patchseq_hash', 'str', 'Hash of patchseq results from amplification and mapping used for updating', {'index': True}),
    # These values are pulled from amplification report
    ('result_BA', 'str', 'Pass/Fail', {'index': True}),
    ('area_400_10000bp', 'float', 'Percentage (0-100) of amplified content in the 400-10,000 bp range which is an indication of intact RNA', {'index': True}),
    ('picogreen_yield', 'float', '(pg/ul)', {'index': True}),
    # These values come from mapping report
    ('cluster_detail', 'str', 'Detailed name of last mapped cluster, for class-level nodes this a descriptive name', {'index': True}),
    ('cluster_label', 'str', 'Label of last mapped cluster, numerical for class-level nodes', {'index': True}),
    ('score', 'float', 'Mapping score from 0-1', {'index': True}),
    ('res_index', 'float', 'Resolution of the last mapped cluster from 0-1 with 1 being the terminal leaf', {'index': True}),
    ('top_leaf', 'str', '', {'index': True}),
    ('top_leaf_score', 'float', 'Confidence of top_leaf mapping (0-1)', {'index': True}),
    ('broad_class_label', 'str', 'Mapped class designation', {'index': True}),
    ('sublass_label', 'str', 'Mapped subclass designation', {'index': True}),
    ('quality_score', 'float', '', {'index': True}),
    ('norm_marker_sum', 'float', '', {'index': True}),
    ('seurat_cluster', 'str', 'Mapped cluster based on Seurat method', {'index': True}),
    ('seurat_score', 'float', 'Mapping score of seurat_cluster (0-1)', {'index': True}),
    # For Tree clustering method, the cummulative scores should = 1
    ('tree_first_cluster', 'str', 'First mapping cluster based on Tree method', {'index': True}),
    ('tree_first_score', 'float', 'Mapping score of first cluster (0-1)', {'index': True}),
    ('tree_second_cluster', 'str', 'Second mapping cluster based on Tree method', {'index': True}),
    ('tree_second_score', 'float', 'Mapping score of second cluster (0-1)', {'index': True}),
    ('tree_third_cluster', 'str', 'Third mapping cluster based on Tree method', {'index': True}),
    ('tree_third_score', 'float', 'Mapping score of third cluster (0-1)', {'index': True}),
    ('tree_call', 'str', 'Tree mapping', {'index': True}),
])

Cell.patch_seq = relationship(PatchSeq, back_populates="cell", cascade="delete", single_parent=True, uselist=False)
PatchSeq.cell = relationship(Cell, back_populates="patch_seq", single_parent=True)