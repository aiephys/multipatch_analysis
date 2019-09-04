from .database import default_db as db

def get_sorted_pulse_responses(pair, query=None):
    """Return all PulseResponse records for *pair* that were generated in a standard induction/recovery
    pulse train, sorted by clamp mode, induction frequency, recovery delay, and pulse number.
    
    Returns
    -------
    sorted_recs : dict
        {('ic', 50.0, 250e-3): {
            1: [pr, pr, pr, ...], 
            2: [pr, pr, ...], 
            ...},
         ...}
         
    Examples
    --------
         
    Get all current clamp pulse responses that are the first pulse in a 50Hz, 250ms delay pulse train:
    
        prs = get_sorted_pulse_responses['ic', 50, 250e-3][1]
    """
    if query is None:
        query = pulse_response_query(pair)

    pr_recs = query.all()

    # group records by (clamp mode, ind_freq, rec_delay), and then by pulse number
    sorted_recs = {}
    for rec in pr_recs:
        stim_key = (rec.patch_clamp_recording.clamp_mode, rec.multi_patch_probe.induction_frequency, rec.multi_patch_probe.recovery_delay)
        sorted_recs.setdefault(stim_key, {})
        pn = rec.stim_pulse.pulse_number
        sorted_recs[stim_key].setdefault(pn, [])
        sorted_recs[stim_key][pn].append(rec)
        
    return sorted_recs


def pulse_response_query(pair, session=None):
    if session is None:
        session = db.session()
    q = session.query(db.PulseResponse, db.PulseResponseFit, db.StimPulse, db.PatchClampRecording, db.MultiPatchProbe, db.Synapse, db.PulseResponse.data)
    q = q.join(db.PulseResponseFit, db.PulseResponse.pulse_response_fit)
    q = q.join(db.StimPulse, db.PulseResponse.stim_pulse)
    q = q.join(db.Recording, db.PulseResponse.recording)
    q = q.join(db.PatchClampRecording, db.PatchClampRecording.recording_id==db.Recording.id)
    q = q.join(db.MultiPatchProbe, db.MultiPatchProbe.patch_clamp_recording_id==db.PatchClampRecording.id)
    q = q.join(db.Synapse, db.Synapse.pair_id==db.PulseResponse.pair_id)
    q = q.filter(db.PulseResponse.pair_id==pair.id)
    return q
