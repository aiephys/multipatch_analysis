from collections import OrderedDict
from pyqtgraph import toposort
from .module import AnalysisModule


def all_modules():
    """Return an ordered dictionary of {module_name:module_class} pairs, sorted by order of dependencies.
    """
    subclasses = AnalysisModule.__subclasses__()
    deps = {c:c.dependencies for c in subclasses}
    return OrderedDict([(mod.name, mod) for mod in toposort(deps)])

