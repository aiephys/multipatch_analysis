from ..pipeline import Pipeline

from .opto_slice import OptoSlicePipelineModule
from .opto_experiment import OptoExperimentPipelineModule
from .opto_morphology import OptoMorphologyPipelineModule
from .opto_cortical_location import OptoCortexLocationPipelineModule
from .opto_dataset import OptoDatasetPipelineModule

class OptoPipeline(Pipeline):
    
    name = 'opto'
    module_classes = [
        OptoSlicePipelineModule,
        OptoExperimentPipelineModule,
        OptoCortexLocationPipelineModule,
        OptoDatasetPipelineModule,
        #OptoMorphologyPipelineModule,
        #PulseResponsePipelineModule,
        #DynamicsPipelineModule,
        #ConnectionStrengthPipelineModule,
        #AverageFirstPulseFitPipelineModule,
        #SingleFirstPulseFitPipelineModule,
    ]
    
    def __init__(self, database, config):
        self.config = config
        self.database = database
        Pipeline.__init__(self, config=config, database=database)