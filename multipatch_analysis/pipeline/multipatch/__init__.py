from ..pipeline import Pipeline

from .slice import SlicePipelineModule
from .experiment import ExperimentPipelineModule
from .dataset import DatasetPipelineModule
from .morphology import MorphologyPipelineModule
from .synapse import SynapsePipelineModule
from .pulse_response import PulseResponsePipelineModule
from .dynamics import DynamicsPipelineModule
from .synapse_prediction import SynapsePredictionPipelineModule
from .resting_state import RestingStatePipelineModule



class MultipatchPipeline(Pipeline):
    
    name = 'multipatch'
    
    module_classes = [
        SlicePipelineModule,
        ExperimentPipelineModule,
        DatasetPipelineModule,
        # MorphologyPipelineModule,
        SynapsePipelineModule,
        PulseResponsePipelineModule,
        SynapsePredictionPipelineModule,
        RestingStatePipelineModule,
        DynamicsPipelineModule,
    ]
    
    def __init__(self, database, config):
        self.config = config
        self.database = database
        Pipeline.__init__(self, config=config, database=database)
