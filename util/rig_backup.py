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
import os, sys, argparse, logging
from aisynphys import util
from aisynphys import config


parser = argparse.ArgumentParser()
parser.add_argument('--jobs', type=str, default="*", help="The name of the backup job(s) to run (default is all jobs listed in config.backup_paths)")
parser.add_argument('--test', action='store_true', default=False, help="Print actions to be taken, do not change any files")
parser.add_argument('--verbose', action='store_true', default=False, help="Verbose output; show files that are skipped over")

args = parser.parse_args(sys.argv[1:])

if args.verbose:
    util.stderr_log_handler.setLevel(logging.DEBUG)

if args.jobs == '*':
    jobs = config.backup_paths.keys()
else:
    jobs = args.jobs.split(',')

for job in jobs:
    spec = config.backup_paths[job]
    source_path = spec['source']
    dest_path = spec['dest']
    log_file = os.path.join(dest_path, 'backup.log')
    util.sync_dir(source_path, dest_path, test=args.test, log_file=log_file)
