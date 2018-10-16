"""
For synchronizing data from archive storage on each rig to backup.

Requires a structure in config.yml like:

backup_paths:
    rig1:
        source: "V:\\"
        dest: "L:\\rig1_backup"
        exclude:
            - "V:\\$Recycle Bin"
    rig2:
        source: "W:\\"
        dest: "L:\\rig2_backup"

"""
from __future__ import print_function
import os, sys, argparse
from multipatch_analysis.util import sync_dir, logger
from multipatch_analysis import config


parser = argparse.ArgumentParser()
parser.add_argument('--rigs', type=str, default="*", help="The name of the rig to back up (default is all rigs)")
parser.add_argument('--test', action='store_true', default=False, help="Print actions to be taken, do not change any files")

args = parser.parse_args(sys.argv[1:])

if args.rigs == '*':
    rigs = config.backup_paths.keys()
else:
    rigs = args.rigs.split(',')

for rig in rigs:
    spec = config.backup_paths[rig]
    source_path = spec['source']
    dest_path = spec['dest']
    log_file = os.path.join(dest_path, 'backup.log')
    sync_dir(source_path, dest_path, test=args.test, log_file=log_file)
