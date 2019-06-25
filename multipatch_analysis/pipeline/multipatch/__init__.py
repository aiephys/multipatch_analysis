from ..pipeline import Pipeline

from .slice import SlicePipelineModule
from .experiment import ExperimentPipelineModule
from .dataset import DatasetPipelineModule
from .morphology import MorphologyPipelineModule
from .pulse_response import PulseResponsePipelineModule
from .dynamics import DynamicsPipelineModule
from .connection_strength import ConnectionStrengthPipelineModule
from .first_pulse_fit import AverageFirstPulseFitPipelineModule
from .first_pulse_fit import SingleFirstPulseFitPipelineModule



class MultipatchPipeline(Pipeline):
    
    name = 'multipatch'
    
    module_classes = [
        SlicePipelineModule,
        ExperimentPipelineModule,
        DatasetPipelineModule,
        MorphologyPipelineModule,
        PulseResponsePipelineModule,
        DynamicsPipelineModule,
        ConnectionStrengthPipelineModule,
        AverageFirstPulseFitPipelineModule,
        SingleFirstPulseFitPipelineModule,
    ]
    
    def __init__(self, database, config):
        self.config = config
        self.database = database
        Pipeline.__init__(self)
        
