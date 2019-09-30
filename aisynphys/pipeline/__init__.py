from .pipeline import Pipeline
from . import multipatch


def all_pipelines():
    """Return a dictionary of {pipeline_name:pipeline_class} pairs
    """
    return Pipeline.all_pipelines()

