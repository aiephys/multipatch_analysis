from neuroanalysis.miesnwb import MiesNwb, MiesSyncRecording, MiesRecording


class MultiPatchExperiment(MiesNwb):
    """Extension of neuroanalysis data abstraction layer to include
    multipatch-specific metadata.
    """
    def create_sync_recording(self, sweep_id):
        return MultiPatchSyncRecording(self, sweep_id)

        
class MultiPatchSyncRecording(MiesSyncRecording):
    def __init__(self, nwb, sweep_id):
        MiesSyncRecording.__init__(self, nwb, sweep_id)
        try:
            self.meta['temperature'] = self.recordings[0].meta['notebook']['Async AD 1: Bath Temperature']
        except Exception:
            self.meta['temperature'] = None
    
    def create_recording(self, sweep_id, ch):
        miesrec = MiesRecording(self, sweep_id, ch)
        stim = miesrec.meta['stim_name'].lower()
        if 'pulsetrain' in stim or 'recovery' in stim:
            return MultiPatchProbe(miesrec)
        else:
            return miesrec


class MultiPatchProbe(MiesRecording):
    def __init__(self, recording):
        self._parent_rec = recording
        
    #@property    
    #def induction_frequency(self):
        #return self.stim_params

    def __getattr__(self, attr):
        return getattr(self._parent_rec, attr)
