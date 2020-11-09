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
from .patch_seq import PatchSeqPipelineModule
from .gap_junction import GapJunctionPipelineModule
from .intrinsic import IntrinsicPipelineModule
from .cortical_location import CortexLocationPipelineModule


class MultipatchPipeline(Pipeline):
    
    name = 'multipatch'
    
    module_classes = [
        SlicePipelineModule,
        ExperimentPipelineModule,
        DatasetPipelineModule,
        IntrinsicPipelineModule,
        PatchSeqPipelineModule,
        MorphologyPipelineModule,
        SynapsePipelineModule,
        GapJunctionPipelineModule,
        PulseResponsePipelineModule,
        #SynapsePredictionPipelineModule,
        RestingStatePipelineModule,
        DynamicsPipelineModule,
        CortexLocationPipelineModule
    ]
    
    def __init__(self, database, config):
        self.config = config
        self.database = database
        Pipeline.__init__(self, config=config, database=database)
