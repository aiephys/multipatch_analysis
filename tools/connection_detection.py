import sys
import user
import pyqtgraph as pg
pg.dbg()

from neuroanalysis.miesnwb import MiesNwb
from multipatch_analysis.experiment_list import cached_experiments
from multipatch_analysis.connection_detection import plot_response_averages

arg = sys.argv[1]
try:
    expt_ind = arg
    all_expts = cached_experiments()
    expt = all_expts[expt_ind].data
except ValueError:
    expt_file = arg
    expt = MiesNwb(expt_file)

plots = plot_response_averages(expt, show_baseline=True, clamp_mode='ic', min_duration=25e-3, pulse_ids=None)

# detect_connections(expt)

if sys.flags.interactive == 0:
    pg.mkQApp().exec_()
