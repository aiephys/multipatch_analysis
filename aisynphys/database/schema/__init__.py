from sqlalchemy.ext.declarative import declarative_base
from ..database import make_table_docstring, make_table as orig_make_table

# schema version should be incremented whenever the schema has changed
schema_version = "15"

# all time series data are downsampled to this rate in the DB
default_sample_rate = 20000
sample_rate_str = '%dkHz' % (default_sample_rate // 1000)


ORMBase = declarative_base()

def make_table(*args, **kwds):
    kwds['ormbase'] = ORMBase
    return orig_make_table(*args, **kwds)

from .meta import *
from .pipeline import *
from .slice import *
from .experiment import *
from .morphology import *
from .dataset import *
from .synapse import *
from .pulse_response import *
from .dynamics import *
from .synapse_prediction import *
from .resting_state_fit import *
from .gap_junction import *

# Create all docstrings now that relationships have been declared
for cls in ORMBase.__subclasses__():
    cls.__doc__ = make_table_docstring(cls)
