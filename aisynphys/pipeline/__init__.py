from .pipeline import Pipeline
from . import pipeline_module
from . import multipatch

from neuroanalysis.util.optional_import import optional_import
opto = optional_import('.opto', package=__name__)


def all_pipelines():
    """Return a dictionary of {pipeline_name:pipeline_class} pairs
    """
    return Pipeline.all_pipelines()

