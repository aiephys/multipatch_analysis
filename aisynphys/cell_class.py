
from __future__ import print_function, division

from sqlalchemy.orm import aliased
from collections import OrderedDict
from .database import default_db
from . import constants


class CellClass(object):
    """Represents a class of cells as a list of selection criteria.
    
    Construct with an arbitrary set of keyword arguments, where each argument specifies
    a criteria for matching cells. Keyword argument names must be a column from the 
    :class:`Cell <aisynphys.database.schema.Cell>`, :class:`Morphology <aisynphys.database.schema.Morphology>`,
    Intrinsic, or CorticalCellLocation database tables.
    Can also be filtered with an arbitrary list of expressions (unnamed arguments),
    each of which must be an sqlalchemy BinaryExpression, as used in query.filter(),
    referring to one of the available tables.
    
    Example::
    
        pv_class = CellClass(cre_type='pvalb')
        inhibitory_class = CellClass(cre_type=('pvalb', 'sst', 'vip'))
        l23_pyr_class = CellClass(pyramidal=True, target_layer='2/3')
        l5_spiny_class = CellClass(dendrite_type='spiny', cortical_layer='5')
        deep_l3_class = CellClass(
            db.CorticalCellLocation.fractional_layer_depth < 0.5,
            cortical_layer='3')
    """
    def __init__(self, name=None, *exprs, **criteria):
        self.criteria = criteria
        self.exprs = exprs
        self._name = name

    @property
    def name(self):
        """A short string representation of this cell class.
        
        If no name was supplied, then this value is `as_tuple` concatenated with spaces.
        """
        if self._name is not None:
            return self._name
        return ' '.join([str(x) for x in self.as_tuple])

    @property
    def as_tuple(self):
        """A tuple representation of this cell class used for display purposes.
        
        Order of elements in the tuple is (target_layer, pyramidal, cre_type), but
        elements are only present if they were specified as criteria for the cell class.
        """
        name = []

        target_layer = self.criteria.get('target_layer')
        cortical_layer = self.criteria.get('cortical_layer')

        if target_layer is not None:
            name.append('L' + target_layer)
        elif cortical_layer is not None:
            name.append('L' + str(cortical_layer))

        if 'dendrite_type' in self.criteria:
            name.append(str(self.criteria['dendrite_type']))

        if 'pyramidal' in self.criteria:
            name.append('pyr' if self.criteria['pyramidal'] else 'nonpyr')

        cre_type = self.criteria.get('cre_type')
        if cre_type is not None:
            name.append(str(cre_type))
        t_type = self.criteria.get('t_type')
        if t_type is not None:
            name.append(str(t_type))
        
        return tuple(name)

    @property
    def is_excitatory(self):
        """True if this is an excitatory cell class, as determined either by
        cre type or morphology.
        """
        cre = self.criteria.get('cre_type')
        if isinstance(cre, tuple):
            cre_is_exc = all([c in constants.EXCITATORY_CRE_TYPES for c in cre])
        else:
            cre_is_exc = cre in constants.EXCITATORY_CRE_TYPES
        pyr = self.criteria.get('pyramidal')
        dendrite = self.criteria.get('dendrite_type')
        return cre_is_exc or pyr is True or dendrite == 'spiny'

    @property
    def output_synapse_type(self):
        """Expected type of synapses "ex" or "in" to be output from this cell type.
        """
        return {True: 'ex', False: 'in'}.get(self.is_excitatory, None)

    def __contains__(self, cell):
        if not (self.criteria or self.exprs):
            return True
        morpho = cell.morphology
        patchseq = cell.patch_seq
        intrinsic = cell.intrinsic
        location = cell.cortical_location
        objs = [cell, morpho, patchseq, intrinsic, location]
        for expr in self.exprs:
            found_attr = False
            key = expr.left.name
            ref_val = expr.right.value
            for obj in objs:
                if hasattr(obj, key):
                    found_attr = True
                    val = getattr(obj, key)
                    if val is not None and expr.operator(val, ref_val):
                        break
                    else:
                        return False
        for k, v in self.criteria.items():
            found_attr = False
            if isinstance(v, dict):
                or_attr = []
                for k2, v2 in v.items():
                    for obj in objs:
                        if hasattr(obj, k2):
                            found_attr = True
                            or_attr.append(getattr(obj, k2) == v2)
                if not any(or_attr):
                    return False
            else:
                for obj in objs:
                    if hasattr(obj, k):
                        found_attr = True
                        if isinstance(v, tuple):
                            if getattr(obj, k) not in v:
                                return False
                        else:
                            if getattr(obj, k) != v:
                                return False
                        break
            if not found_attr:
                # return False
                raise Exception('Cannot use "%s" for cell typing; attribute not found on cell or linked objects' % k)
        return True

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, a):
        """Cell class is considered equal to its *name* to allow it to be indexed from a dict more
        easily::

            cc = CellClass(cre_type='sst', layer='6')
            cc.name => 'L6 sst'
            {cc: 1}['L6 sst'] => 1 
        """
        if isinstance(a, str):
            return a == self.name
        elif isinstance(a, CellClass):
            return a.name == self.name
        else:
            return object.__eq__(self)  # should raise NotImplemented

    def __repr__(self):
        return "<CellClass %s>" % self.name

    def __str__(self):
        return self.name
        
    def filter_query(self, query, cell_table, db=None):
        """Return a modified query (sqlalchemy) that filters results to include only those in
        this cell class.
        """
        if db is None:
            db = default_db
        morpho = aliased(db.Morphology)
        intrinsic = aliased(db.Intrinsic)
        location = aliased(db.CorticalCellLocation)
        query = (query.outerjoin(morpho, morpho.cell_id==cell_table.id)
                 .outerjoin(intrinsic, intrinsic.cell_id==cell_table.id)
                 .outerjoin(location, location.cell_id==cell_table.id))
        tables = [cell_table, morpho, intrinsic, location]
        for expr in self.exprs:
            found_attr = False
            key = expr.left.name
            ref_val = expr.right.value
            for table in tables:
                if hasattr(table, key):
                    found_attr = True
                    query = query.filter(expr.operator(getattr(table, key), ref_val))
                    break
            if not found_attr:
                raise Exception('Cannot use "%s" for cell typing; attribute not found in available tables.' % key)
        for k, v in self.criteria.items():
            found_attr = False
            for table in tables:
                if hasattr(table, k):
                    found_attr = True
                    if isinstance(v, tuple):
                        query = query.filter(getattr(table, k).in_(v))
                    else:
                        query = query.filter(getattr(table, k)==v)
                    break
            if not found_attr:
                raise Exception('Cannot use "%s" for cell typing; attribute not found in available tables.' % k)

        return query                
                
            


def classify_cells(cell_classes, cells=None, pairs=None, missing_attr='raise'):
    """Given cell class definitions and a list of cells, return a dict indicating which cells
    are members of each class.

    Parameters
    ----------
    cell_classes : dict
        List of CellClass instances
    cells : list | None
        List of Cell instances to be classified.
    pairs : list | None
        List of pairs from which cells will be collected. May not be used with *cells*
    missing_attr : str
        Determines the behavior when a criteria attribute is missing on a 
        cell. If 'ignore', then the cell is excluded from the result,. If 'raise',
        then an exception is raised. Default is 'raise'.
        
    Returns
    -------
    cell_groups : OrderedDict
        Dictionary mapping {cell_class: [list of cells]}
        
        
    Example
    -------
        
        pv_cell_class = CellClass(cre_type='pvalb', target_layer='2/3')
        sst_cell_class = CellClass(cre_type='sst', target_layer='2/3')
        
        cell_classes = [pv_cell_class, sst_cell_class]
        cells = session.Query(db.Cell).all()
        
        grouped_cells = classify_cells(cell_classes, cells=cells)
    
        pv_cells = grouped_cells[pv_cell_class]    
        sst_cells = grouped_cells[sst_cell_class]    

    """
    if pairs is not None:
        assert cells is None, "cells and pairs arguments are mutually exclusive"
        cells = set([p.pre_cell for p in pairs] + [p.post_cell for p in pairs])
    cell_groups = OrderedDict([(cell_class, set()) for cell_class in cell_classes])
    for cell in cells:
        for cell_class in cell_classes:
            try:
                if cell in cell_class:
                    cell_groups[cell_class].add(cell)
            except Exception:
                if missing_attr == 'ignore':
                    continue
                else:
                    raise
    return cell_groups


def classify_pairs(pairs, cell_groups):
    """Given a list of cell pairs and a dict that groups cells together by class (ie the output of classify_cells),
    return a dict that groups pairs into (pre, post) cell type buckets.
    
    Parameters
    ----------
    pairs : list of Pair instances
        The Pair instances (probably returned from a database query) to be grouped
    cell_groups : dict
        Specifies the cell classes and the cells that belong to each class. The format is
        the same as the output of classify_cells().
        
    Returns
    -------
    pair_groups : OrderedDict
        Maps {(pre_class, post_class): [list of pairs]}
    """
    results = OrderedDict()
    for pre_class, pre_group in cell_groups.items():
        for post_class, post_group in cell_groups.items():
            post_group = cell_groups[post_class]
            class_pairs = [p for p in pairs if p.pre_cell in pre_group and p.post_cell in post_group]
            results[(pre_class, post_class)] = class_pairs
    
    return results
