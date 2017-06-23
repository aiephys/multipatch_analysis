from neuroanalysis.miesnwb import MiesNwb, MiesSyncRecording, MiesRecording


class MultipatchExperiment(MiesNwb):
    """Extension of neuroanalysis data abstraction layer to include
    multipatch-specific metadata.
    """
    def create_sync_recording(self, sweep_id):
        return MultipatchSyncRecording(self, sweep_id)

        
class MultipatchSyncRecording(MiesSyncRecording):
    def __init__(self, nwb, sweep_id):
        MiesSyncRecording.__init__(self, nwb, sweep_id)
        self.meta['temperature'] = self.recordings[0].meta['notebook']['Async AD 1: Bath Temperature']
    
    def create_recording(self, sweep_id, ch):
        return MiesRecording(self, sweep_id, ch)



