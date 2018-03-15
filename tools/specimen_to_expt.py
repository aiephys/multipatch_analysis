from multipatch_analysis.experiment_list import cached_experiments
import argparse
import sys
import pyqtgraph as pg

specimen_id = sys.argv[1]

expts = cached_experiments()
expt = [e for e in expts if e.specimen_id == specimen_id]
for site in expt:
    print ('For specimen %s:\t \n Backup path:\t %s \n Original path:\t %s \n Rig:\t %s \n UID:\t %s' % (site.specimen_id, site.path, site.original_path, site.rig_name, site.uid))