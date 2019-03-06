from .analysis_module import AnalysisModule
from .slice import SliceAnalysisModule


def all_modules():
    """Return an ordered dictionary of {module_name:module_class} pairs, sorted by order of dependencies.
    """
    return AnalysisModule.all_modules()

