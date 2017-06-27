# *-* coding: utf-8 *-*
"""
CLI script for generating reports from multipatch experiment data.

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


parser = argparse.ArgumentParser()
parser.add_argument('--region', type=str)
parser.add_argument('--start', type=arg_to_date)
parser.add_argument('--stop', type=arg_to_date)
parser.add_argument('--list-stims', action='store_true', default=False, dest='list_stims',
                    help='print a list of each connection and the stimulus sets acquired')
parser.add_argument('--sweep-threshold', nargs = '*', type=int, action='store', default=[5,10], dest='sweep_threshold',
                    help='Combined with --list-stims, for each connection type, prints the number of connections'
                         '' 'for which there are >= sweep_threshold number of sweeps/stimulus set. Two thresholds'
                         '' 'are set one for induction protocols (default=5) and one for recovery (default=10')
parser.add_argument('files', nargs='*', type=os.path.abspath)
parser.add_argument('--cre_type', nargs=2, type=str)
parser.add_argument('--calcium', type=str, help='define external calcium concentration as "Low" or "High"')
parser.add_argument('--age', nargs=2, type=int, help='Define age as a range from min to max')
parser.add_argument('--temp', type=int)
args = parser.parse_args(sys.argv[1:])

cache_file = 'expts_cache.pkl'
all_expts = ExperimentList(cache=cache_file)

for f in args.files:
    all_expts.load(f)

if len(all_expts) == 0:
    print("No experiments loaded; bailing out.")
    sys.exit(-1)

expts = all_expts.select(start=args.start, stop=args.stop, region=args.region, cre_type=args.cre_type, calcium=args.calcium, age=args.age, temp=args.temp)
if len(args.files) > 0:
    expts = expts.select(source_files=args.files)

if len(expts) == 0:
    print("No experiments selected; bailing out.")
    sys.exit(-1)

expts.check()


# Print list of experiments
expts.print_expt_summary(args.list_stims)

# Print list of connections found
expts.print_connection_summary(args.cre_type, args.list_stims)

# Print stimulus summary for each connection type
if args.list_stims:
    expts.print_connection_sweep_summary(args.sweep_threshold)

# Generate a summary of connectivity
expts.print_connectivity_summary(args.cre_type)

# Print extra information about labeling
expts.print_label_summary()


pg.mkQApp()

plots = expts.distance_plot('sim1', 'sim1', color=(0, 150, 255))
expts.distance_plot('tlx3', 'tlx3', plots=plots, color=(200, 100, 0))
#expts.distance_plot('pvalb', 'pvalb', plot=p, color=(200, 0, 200))

types = ['unknown', 'sim1', 'tlx3', 'pvalb', 'sst', 'vip']
#types = ['sim1', 'unknown']
expts.matrix(types, types)

# cache everything!
all_expts.write_cache()
