from .pipeline import Pipeline
from . import pipeline_module
from . import multipatch
from . import opto


def all_pipelines():
    """Return a dictionary of {pipeline_name:pipeline_class} pairs
    """
    return Pipeline.all_pipelines()

