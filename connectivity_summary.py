# *-* coding: utf-8 *-*
"""
I have a summary of ALM multipatch experiments stored on workflowy and exported to expt_summary.txt.

This file describes, for each experiment:
   - which cre labels were used
   - which cells are +/- for which labels
   - which cells had usable recordings
   - which cells were connected

The purpose of this script is to parse the exported file and return a summary
of total connections detected / probed for each cell type pair.

"""
from __future__ import print_function, division

import argparse
import datetime
import os
import re
import sys

import pyqtgraph as pg

from experiment_list import ExperimentList


def arg_to_date(arg):
    if arg is None:
        return None
    parts = re.split('\D+', arg)
    return datetime.date(*map(int, parts))




# read exported file

#steph = open('/home/luke/mnt/mp1/Steph/DataSummary', 'r').readlines()
#pasha = open('/home/luke/mnt/mp2/data/Pasha/connection_analysis', 'r').readlines()
#alex = open('/home/luke/mnt/mp3/data/Alex/connection analysis', 'r').readlines()
#lines = steph + pasha + alex



# first parse indentation to generate a hierarchy of strings
# Note: summary was exported in plain text mode, which can be difficult to parse.
# Might want to switch to XML if this becomes an issue.



# Now go through each experiment and read cell type / connectivity data


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--region', type=str)
    parser.add_argument('--start', type=arg_to_date)
    parser.add_argument('--stop', type=arg_to_date)
    parser.add_argument('--list-stims', action='store_true', default=False, dest='list_stims',
                        help='print a list of each connection and the stimulus sets acquired')
    parser.add_argument('--sweep-threshold', type=int, action='store', default=5, dest='sweep_threshold',
                        help='Combined with --list-stims, for each connection type, prints the number of connections'
                             '' 'for which there are >= sweep_threshold number of each stimulus set')
    parser.add_argument('files', nargs='*', type=os.path.abspath)
    args = parser.parse_args(sys.argv[1:])

    cache_file = 'expts_cache.pkl'
    all_expts = ExperimentList(cache=cache_file)

    for f in args.files:
        all_expts.load(f)

    if len(all_expts) == 0:
        print("No experiments loaded; bailing out.")
        sys.exit(-1)

    expts = all_expts.select(start=args.start, stop=args.stop, region=args.region)
    if len(args.files) > 0:
        expts = expts.select(source_files=args.files)

    if len(expts) == 0:
        print("No experiments selected; bailing out.")
        sys.exit(-1)

    expts.check()
    

    # Print list of experiments
    expts.print_expt_summary(args.list_stims)

    # Print list of connections found
    expts.print_connection_summary(args.list_stims)

    # Print stimulus summary for each connection type
    if args.list_stims:
        expts.print_connection_sweep_summary(args.sweep_threshold)

    # Generate a summary of connectivity
    expts.print_connectivity_summary()

    # Print extra information about labeling
    expts.print_label_summary()


    p = pg.plot()
    expts.distance_plot('sim1', 'sim1', plot=p, color=(0, 0, 255))
    expts.distance_plot('tlx3', 'tlx3', plot=p, color=(200, 200, 0))
    expts.distance_plot('pvalb', 'pvalb', plot=p, color=(200, 0, 200))

    types = ['unknown', 'sim1', 'tlx3', 'pvalb', 'sst', 'vip']
    #types = ['sim1', 'unknown']
    expts.matrix(types, types)
    
    # cache everything!
    all_expts.write_cache()
    