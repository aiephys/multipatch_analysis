from sqlalchemy.ext.declarative import declarative_base
from ..database import make_table_docstring, make_table as orig_make_table

# schema version should be incremented whenever the schema has changed
schema_version = "19"

# all time series data are downsampled to this rate in the DB
default_sample_rate = 20000
sample_rate_str = '%dkHz' % (default_sample_rate // 1000)


ORMBase = declarative_base()


def schema_description():
    """Return a structure describing tables and columns in the schema.
    """
    schema = {}

    for table in ORMBase.metadata.sorted_tables:
        cols = {}
        for colname,col in table.columns.items():
            try:
                pytype = col.type.python_type
            except Exception:
                pytype = None
            cols[colname] = {'type': str(col.type), 'pytype': pytype, 'comment': col.comment}
        schema[table.name] = {'columns': cols, 'comment': table.comment}

    return schema


def make_table(*args, **kwds):
    kwds['ormbase'] = ORMBase
    return orig_make_table(*args, **kwds)

from .meta import *
from .pipeline import *
from .slice import *
from .experiment import *
from .intrinsic import *
from .morphology import *
from .dataset import *
from .synapse import *
from .pulse_response import *
from .dynamics import *
from .synapse_prediction import *
from .resting_state_fit import *
from .gap_junction import *
from .cortical_location import *
from .patch_seq import *
from .synapse_model import *

# Create all docstrings now that relationships have been declared
for cls in ORMBase.__subclasses__():
    cls.__doc__ = make_table_docstring(cls)
