from __future__ import division
import time, datetime
import aisynphys.database.database as db
from neuroanalysis.ui.plot_grid import PlotGrid
s = db.session()


q = """
    select 
        pcrec.baseline_rms_noise,
        substring(experiment.original_path from 36 for 1),
        recording.device_key,
        recording.start_time
    from 
        patch_clamp_recording pcrec
        join recording on pcrec.recording_id=recording.id
        join sync_rec on recording.sync_rec_id=sync_rec.id
        join experiment on sync_rec.experiment_id=experiment.id
    where
        pcrec.clamp_mode='ic'
        and pcrec.baseline_rms_noise is not null
        and recording.device_key is not null
        and experiment.original_path is not null
"""


rec = s.execute(q)
rows = rec.fetchall()

import pyqtgraph as pg
import numpy as np

rms = np.array([row[0] for row in rows])
rig = np.array([row[1] for row in rows]).astype(int)
hs = np.array([row[2] for row in rows]).astype(int)
col = rig*8 + hs

ts = np.array([time.mktime(row[3].timetuple()) for row in rows])
#ts -= ts[0]

pg.plot(col + np.random.uniform(size=len(col))*0.7, rms, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 255, 255, 50))

plt = pg.plot(labels={'left': 'number of sweeps (normalized per rig)', 'bottom': ('baseline rms error', 'V')})
plt.addLegend()

grid = PlotGrid()
grid.set_shape(3, 1)
grid.show()

for r, c in ((1, 'r'), (2, 'g'), (3, 'b')):
    mask = rig==r
    rig_data = rms[mask]
    y, x = np.histogram(rig_data, bins=np.linspace(0, 0.002, 1000))
    plt.plot(x, y/len(rig_data), stepMode=True, connect='finite', pen=c, name="Rig %d" % r)

    p = grid[r-1, 0]
    p.plot(ts[mask], rig_data, pen=None, symbol='o', symbolPen=None, symbolBrush=(255, 255, 255, 100))
    p.setLabels(left=('rig %d baseline rms noise'%r, 'V'))

grid.setXLink(grid[0, 0])
grid.setYLink(grid[0, 0])

