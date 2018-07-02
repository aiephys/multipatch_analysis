from __future__ import print_function
import io
import numpy as np
import pyqtgraph as pg
from neuroanalysis.data import Trace, TraceList
import datetime
import multipatch_analysis.database.database as db


session = db.Session()

ts = datetime.datetime.fromtimestamp(1502312765.005)
pre_ch = 0
post_ch = 6

q = session.query(db.Experiment.id).filter(db.Experiment.acq_timestamp==ts)
exp = q[0]


conditions = [
    "experiment.id=%s" % exp.id,
    "pre_rec.device_key=%s" % pre_ch,
    "post_rec.device_key=%s" % post_ch,
    "post_pcrec.clamp_mode='ic'",
    "multi_patch_probe.induction_frequency<100",
    "(stim_pulse.pulse_number=1 or stim_pulse.pulse_number=9)",
]


query = """
    select 
        pulse_response.data 
    from 
        pulse_response
        join stim_pulse on pulse_response.stim_pulse_id=stim_pulse.id
        join recording post_rec on pulse_response.recording_id=post_rec.id
        join patch_clamp_recording post_pcrec on post_pcrec.recording_id=post_rec.id
        join multi_patch_probe on multi_patch_probe.patch_clamp_recording_id=post_pcrec.id
        join recording pre_rec on stim_pulse.recording_id=pre_rec.id
        join sync_rec on post_rec.sync_rec_id=sync_rec.id
        join experiment on sync_rec.experiment_id=experiment.id
    where
        {conditions}
""".format(conditions='\n        and '.join(conditions))

print(query)

rp = session.execute(query)

recs = rp.fetchall()
data = [np.load(io.BytesIO(rec[0])) for rec in recs]
print("\n\nloaded %d records" % len(data))



plt = pg.plot(labels={'left': ('Vm', 'V')})
traces = TraceList()
for i,x in enumerate(data):
    trace = Trace(x - np.median(x[:100]), sample_rate=20000)
    traces.append(trace)
    if i<100:
        plt.plot(trace.time_values, trace.data, pen=(255, 255, 255, 100))

avg = traces.mean()
plt.plot(avg.time_values, avg.data, pen='g')

