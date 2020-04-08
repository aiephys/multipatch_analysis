### Utility functions for extracting recordings from the NWB
import numpy as np

def get_lp_sweeps(sweeps, dev_id):
    # get stimulus sweeps from the TargetV and IF Curve stim sets
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
    # returns the start and stop time of the target_v or If_curve long pulse
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
    trodes = {e.device_id: e.id for e in expt.electrodes}
    trode_id = trodes[recording.device_id]
    srecs = {srec.ext_id: srec for srec in expt.sync_recs}
    try:
        srec = srecs[recording.parent.key]
    except KeyError:
        return None
    recs = {rec.electrode_id:rec for rec in srec.recordings}
    return recs.get(trode_id, None)