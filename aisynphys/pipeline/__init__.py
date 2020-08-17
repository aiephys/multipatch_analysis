from .pipeline import Pipeline
from . import pipeline_module
from . import multipatch
try:
    from . import opto
except ImportError:
    from aisynphys.util import optional_import
    opto = optional_import('aisynphys.pipeline.opto')


def all_pipelines():
    """Return a dictionary of {pipeline_name:pipeline_class} pairs
    """
    return Pipeline.all_pipelines()

