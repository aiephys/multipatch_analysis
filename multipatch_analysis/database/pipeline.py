from collections import OrderedDict
from .database import TableGroup

__all__ = ['pipeline_tables', 'Pipeline']


class PipelineTableGroup(TableGroup):
    """Tables used for storing metadata about pipeline job status.
    """
    schemas = OrderedDict([
        ('pipeline', [
            "Stores information about which pipeline analysis jobs were run, when, and whether there was an error.",
            ('module_name', 'str', 'The name of the pipeline module that generated this result', {'index': True}),
            ('job_id', 'float', 'Unique value identifying the job that was processed', {'index': True}),
            ('finish_time', 'datetime', 'The date/time when this job completed processing'),
            ('success', 'bool', 'Whether the job completed successfully', {'index': True}),
            ('error', 'str', 'Error or warning messages generated during job processing'),
        ]),
    ])

pipeline_tables = PipelineTableGroup()
Pipeline = pipeline_tables['pipeline']
