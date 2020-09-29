### Utility functions for extracting recordings from the NWB
import numpy as np
from collections import defaultdict

def get_intrinsic_recording_dict(expt, dev_id, check_qc=True):
    """Get stimulus recordings for intrinsic stimuli
    """
    recording_dict = defaultdict(list)
    sweeps = expt.data.contents
    for sweep in sweeps:
        devices = sweep.devices
        if dev_id not in devices:
            continue
        rec = sweep[dev_id]
        stim_name = rec.stimulus.description
        sweep_types = ['TargetV', 'If_Curve', 'IV_Curve ', 'Chirp']
        if rec.clamp_mode=='ic' and any([code in stim_name for code in sweep_types]):
            if check_qc:
                db_rec = get_db_recording(expt, rec)
                if db_rec is None or db_rec.patch_clamp_recording.qc_pass is False:
                    continue
            code = "Chirp" if 'Chirp' in stim_name else "LP" 
            recording_dict[code].append(rec)
    return recording_dict

def get_lp_sweeps(sweeps, dev_id):
    """Get stimulus sweeps from the TargetV (subthreshold, mostly hyperpolarizing) 
    and IF Curve (suprathreshold) stim sets
    """
    targv_sweeps = []
    ifc_sweeps = []
    for sweep in sweeps:
        devices = sweep.devices
        if dev_id not in devices:
            continue
        rec = sweep[dev_id]
        stim_name = rec.stimulus.description
        targv = 'TargetV' in stim_name
        ifc = 'If_Curve' in stim_name
        if targv:
            targv_sweeps.append(sweep)
        if ifc:
            ifc_sweeps.append(sweep)
    return targv_sweeps, ifc_sweeps

def get_pulse_times(recording):
    """Returns the start and stop time of the target_v or If_curve long pulse.
    How this gets parsed depends on how the experiment was aquired as these stimuli 
    were added later.
    """
    stims = recording.stimulus.items
    # newer experiments have 5 elements, the 4th one is the long pulse of the target_v or If_curve stimulus
    if len(stims) == 5:
        return stims[3].start_time, stims[3].start_time + stims[3].duration
    # older experiments only have metadata for the test pulse, use the following to parse the 
    # wave itself to find the stimlus
    elif len(stims) == 2:
        # skip over test pulse
        start = stims[1].start_time + stims[1].duration
        cmd = recording['command'].time_slice(start, None)
        # find pulse by parsing command wave
        dif = np.diff((cmd.data != cmd.data[0]).astype(int))
        inds = np.argwhere(dif != 0)[:, 0]
        if len(inds) == 0:
            return None
        return cmd.time_at(inds[0]), cmd.time_at(inds[1])
    else:
        return None

def get_db_recording(expt, recording):
    """Get the database record for the recording given the raw recording object
    """
    trodes = {e.device_id: e.id for e in expt.electrodes}
    trode_id = trodes[recording.device_id]
    srecs = {srec.ext_id: srec for srec in expt.sync_recs}
    try:
        srec = srecs[recording.parent.key]
    except KeyError:
        return None
    recs = {rec.electrode_id:rec for rec in srec.recordings}
    return recs.get(trode_id, None)