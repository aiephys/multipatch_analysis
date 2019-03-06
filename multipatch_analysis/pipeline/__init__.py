from .pipeline_module import PipelineModule
from .slice import SlicePipelineModule


def all_modules():
    """Return an ordered dictionary of {module_name:module_class} pairs, sorted by order of dependencies.
    """
    return PipelineModule.all_modules()

