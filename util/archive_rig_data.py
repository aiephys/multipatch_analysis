"""
Used to synchronize data on each rig from the primary storage drive (fast/small) to the archive storage drive (slow/large).
"""
import sys, argparse
from multipatch_analysis import config, util


def archive_rig_data(rig_name, test=False):
    paths = config.rig_data_paths[rig_name]
    for path_set in paths:
        source = path_set['primary']
        dest = path_set['archive']
        print("Archiving data for %s : %s => %s..." % (rig_name, source, dest))
        util.sync_dir(source, dest, test=test)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rig', type=str, default=config.rig_name, help="The name of the rig to archive (must be listed in config.rig_data_paths)")
    parser.add_argument('--test', action='store_true', default=False, help="Print actions to be taken, do not change any files")
    
    args = parser.parse_args(sys.argv[1:])

    if not args.rig[-1].isdigit():
        raise Exception('Rig name "%s" is invalid.')

    archive_rig_data(rig_name=args.rig, test=args.test)
