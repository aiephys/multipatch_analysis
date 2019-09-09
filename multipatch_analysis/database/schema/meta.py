from collections import OrderedDict
from . import make_table

__all__ = ['Metadata']


Metadata = make_table(
    name='metadata',
    comment="A single metadata related to the entire database.",
    columns=[
    ]
)
