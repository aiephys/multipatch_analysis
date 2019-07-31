from collections import OrderedDict
from . import make_table

__all__ = ['Pipeline']


Pipeline = make_table(
    name='pipeline',
    comment="Stores information about which pipeline analysis jobs were run, when, and whether there was an error.",
    columns=[
        ('module_name', 'str', 'The name of the pipeline module that generated this result', {'index': True}),
        ('job_id', 'str', 'Unique value identifying the job that was processed', {'index': True}),
        ('finish_time', 'datetime', 'The date/time when this job completed processing'),
        ('success', 'bool', 'Whether the job completed successfully', {'index': True}),
        ('error', 'str', 'Error or warning messages generated during job processing'),
    ]
)
