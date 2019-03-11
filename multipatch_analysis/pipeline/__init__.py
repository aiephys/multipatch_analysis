from .pipeline_module import PipelineModule
from .slice import SlicePipelineModule
from .experiment import ExperimentPipelineModule
from .dataset import DatasetPipelineModule
from .morphology import MorphologyPipelineModule
from .pulse_response import PulseResponsePipelineModule
from .connection_strength import ConnectionStrengthPipelineModule
from .first_pulse_fit import FirstPulseFitPipelineModule


def all_modules():
    """Return an ordered dictionary of {module_name:module_class} pairs, sorted by order of dependencies.
    """
    return PipelineModule.all_modules()

