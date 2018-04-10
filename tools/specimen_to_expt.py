from multipatch_analysis.experiment_list import cached_experiments
import argparse
import sys
import pyqtgraph as pg

parser = argparse.ArgumentParser()
parser.add_argument('--specimen-id', nargs='*')
parser.add_argument('--uid', nargs='*')

args = parser.parse_args(sys.argv[1:])

expts = cached_experiments()
if args.specimen_id is not None:
    specimen_list = args.specimen_id
    for specimen in specimen_list: 
        expt = [e for e in expts if specimen in e.specimen_id]
        for site in expt:
            print ('For specimen %s:\t \n Backup path:\t %s \n Original path:\t %s \n Rig:\t %s \n UID:\t %s' % (site.specimen_id, site.path, site.original_path, site.rig_name, site.uid))
elif args.uid is not None:
    specimen_list = args.uid
    for specimen in specimen_list:
        expt = expts[args.uid]
        print (' Specimen ID:\t %s \n Backup path:\t %s \n Original path:\t %s \n Rig:\t %s \n UID:\t %s' % (expt.specimen_id, expt.path, expt.original_path, expt.rig_name, expt.uid))