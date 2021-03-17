
from __future__ import print_function, division

from sqlalchemy.orm import aliased
import sqlalchemy.sql.elements
from collections import OrderedDict
from .database import default_db
from .database.schema import schema_description
from . import constants


# tables containing data used to classify cells, mapped to the cell attributes
# used to access the table
_cell_data_tables = [
    ('cell', None), 
    ('morphology', 'morphology'),
    ('patch_seq', 'patch_seq'), 
    ('intrinsic', 'intrinsic'),
    ('cortical_cell_location', 'cortical_location'),
]

# names of attributes available for classification, mapped back to their source tables
_db_schema = schema_description()
_criteria_attributes = {}
for table_name,table_attr in _cell_data_tables:
    cols = _db_schema[table_name]['columns']
    for k in cols.keys():
        _criteria_attributes[k] = (table_name, table_attr)


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
        l23_pyr_class = CellClass(cortical_layer='2/3')
        l5_spiny_class = CellClass(dendrite_type='spiny', cortical_layer='5')
        deep_l3_class = CellClass(
            db.CorticalCellLocation.fractional_layer_depth < 0.5,
            cortical_layer='3')
    """

    def __init__(self, *exprs, name=None, **criteria):        
        # sanity check inputs
        global _criteria_attributes
        assert name is None or isinstance(name, str), f"name must be a string or None (got {repr(name)})"
        for k,v in criteria.items():
            assert k in _criteria_attributes, f"Key '{k}' is not a valid cell class criterion."
        for ex in exprs:
            assert isinstance(ex, sqlalchemy.sql.elements.BinaryExpression), f"non-keyword arguments must be sqlalchemy binary expressions (got {repr(ex)})"
            assert ex.left.name in _criteria_attributes, f"Key '{ex.left.name}' is not a valid cell class criterion."

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
        """True if this class includes only excitatory cells; False if this class includes only inhibitory cells;
        None if the class may include a mixture of excitatory and inhibitory cells.

        Relevant criteria used here are:
        * cell.cre_type
        * morphology.dendrite_type
        * cell.cell_class
        * cell.cell_class_nonsynaptic

        """
        is_ex = []

        cre = self.criteria.get('cre_type')
        if not isinstance(cre, (tuple, list)):
            cre = (cre,)
        cre_is_exc = all([c in constants.EXCITATORY_CRE_TYPES for c in cre])
        cre_is_inh = all([c in constants.INHIBITORY_CRE_TYPES for c in cre])
        if cre_is_exc:
            is_ex.append(True)
        elif cre_is_inh:
            is_ex.append(False)

        dendrite = self.criteria.get('dendrite_type')
        if dendrite == 'spiny':
            is_ex.append(True)
        elif dendrite in ['aspiny', 'sparsely spiny']:
            is_ex.append(False)
        elif dendrite is not None:
            return None

        for cell_class in [self.criteria.get('cell_class'), self.criteria.get('cell_class_nonsynaptic')]:
            if cell_class == 'ex':
                is_ex.append(True)
            elif cell_class == 'in':
                is_ex.append(False)
            elif cell_class == 'mixed':
                return None
            elif cell_class is not None:
                raise ValueError("cell class criteria must be 'ex', 'in', or 'mixed'")

        if len(is_ex) == 0:
            return None

        if all(is_ex):
            return True
        elif not any(is_ex):
            return False
        else:
            return None

    @property
    def output_synapse_type(self):
        """Expected type of synapses "ex", "in", or None to be output from this cell type.
        """
        return {True: 'ex', False: 'in'}.get(self.is_excitatory, None)

    def __contains__(self, cell):
        if len(self.criteria) == 0 and len(self.exprs) == 0:
            return True

        # check expressions
        for expr in self.exprs:
            # get variable name and value from expression
            key = expr.left.name
            ref_val = expr.right.value

            # check requested value
            val = self._get_cell_subattr(cell, key)
            if val is not None and expr.operator(val, ref_val):
                continue
            else:
                return False

        # check keyword arg criteria
        for k, v in self.criteria.items():
            if isinstance(v, dict):
                or_attr = []
                for k2, v2 in v.items():
                    v1 = self._get_cell_subattr(cell, k2)
                    or_attr.append(v1 == v2)
                if not any(or_attr):
                    return False
            elif isinstance(v, (tuple, list)):
                if self._get_cell_subattr(cell, k) not in v:
                    return False
            else:
                if self._get_cell_subattr(cell, k) != v:
                    return False

        return True

    @staticmethod
    def _get_cell_subattr(cell, attr):
        """Get attribute from cell or one of its linked tables (morphology, intrinsic, etc..)
        """
        global _criteria_attributes
        if attr not in _criteria_attributes:
            raise Exception(f'Cannot use "{attr}" for cell typing; attribute not found on cell or linked objects')
        sub_obj_name = _criteria_attributes[attr][1]
        obj = cell if sub_obj_name is None else getattr(cell, sub_obj_name, None)
        if obj is None:
            return None
        return getattr(obj, attr)

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
                    if isinstance(v, (tuple, list)):
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


_criteria_attribute_cache = {}
def _get_criteria_attributes(db):
    """Return a dict mapping attribute:(table_attribute, db_table) for all cell-related attributes
    that can be used for classification criteria.
    """
    global _criteria_attribute_cache

    if db not in _criteria_attribute_cache:
        # tables containing extra per-cell data that can be used as classification criteria
        criteria_tables = [
            (None, db.Cell), 
            ('morphology', db.Morphology), 
            ('patch_seq', db.PatchSeq), 
            ('intrinsic', db.Intrinsic), 
            ('cortical_location', db.CorticalCellLocation),
        ]

        # names of attributes available for classification, mapped back to their source tables
        criteria_attributes = {}
        for name,table in criteria_tables:
            for k in table.__table__.columns.keys():
                criteria_attributes[k] = name,table

        _criteria_attribute_cache[db] = criteria_attributes

    return _criteria_attribute_cache[db]
