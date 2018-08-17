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
from collections import OrderedDict
import numpy as np
import pyqtgraph as pg

from multipatch_analysis.experiment_list import ExperimentList, cache_file
from multipatch_analysis import config

def arg_to_date(arg):
    if arg is None:
        return None
    parts = re.split('\D+', arg)
    return datetime.date(*map(int, parts))


parser = argparse.ArgumentParser()
parser.add_argument('--reload', action='store_true', default=False, dest='reload',
                    help='Reload all experiment data from the server.')
parser.add_argument('--reload-old', action='store_true', default=False, dest='reload_old',
                    help='Reload all experiment data from old summary files.')
parser.add_argument('--region', type=str)
parser.add_argument('--organism', type=str, help='"mouse" or "human"')
parser.add_argument('--start', type=arg_to_date)
parser.add_argument('--stop', type=arg_to_date)
parser.add_argument('--list-stims', action='store_true', default=False, dest='list_stims',
                    help='print a list of each connection and the stimulus sets acquired')
parser.add_argument('--sweep-threshold', nargs = '*', type=int, action='store', default=[5,10], dest='sweep_threshold',
                    help='Combined with --list-stims, for each connection type, prints the number of connections '
                         'for which there are >= sweep_threshold number of sweeps/stimulus set. Two thresholds '
                         'are set one for induction protocols (default=5) and one for recovery (default=10')
parser.add_argument('files', nargs='*', type=os.path.abspath)
parser.add_argument('--cre-type', nargs=2, type=str)
parser.add_argument('--calcium', type=str,
                    help='define external calcium concentration as "Low", "High"')
parser.add_argument('--age', type=str, help='Define age as a range from min to max.  Ex age=30-40')
parser.add_argument('--temp', type=int)

args = parser.parse_args(sys.argv[1:])

all_expts = ExperimentList(cache=cache_file)

if args.reload:
    all_expts.load_from_server()

if args.reload_old:
    files = config.summary_files
    if len(files) == 0:
        print("No old-format summary files given in config.yml:summary_files")
        sys.exit(-1)
    for f in files:
        all_expts.load(f)
    
for f in args.files:
    all_expts.load(f)

if len(all_expts) == 0:
    print("No experiments loaded; bailing out.")
    sys.exit(-1)

# force cell QC to run before caching
for expt in all_expts:
    expt.connections_probed

# cache everything!
all_expts.write_cache()
print("Cache successfully updated!")


for i, ex in enumerate(all_expts):
    ex.summary_id = i

expts = all_expts.select(start=args.start, stop=args.stop, region=args.region, cre_type=args.cre_type, calcium=args.calcium, age=args.age, temp=args.temp, organism=args.organism)
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
    expts.print_connection_sweep_summary(args.cre_type, args.sweep_threshold)

# Generate a summary of connectivity
expts.print_connectivity_summary(args.cre_type)

# Print extra information about labeling
expts.print_label_summary()


pg.mkQApp()


mouse_types = [
    ('1', 'ndnf'),
    ('2/3', 'unknown'),
    ('2/3', 'pvalb'),
    ('2/3', 'sst'),
    ('2/3', 'vip'),
    # ('4', 'unknown'),
    ('4', 'nr5a1'),
    ('4', 'pvalb'),
    ('4', 'sst'),
    ('4', 'vip'),
    # ('5', 'unknown'),
    ('5', 'sim1'),
    ('5', 'tlx3'),
    ('5', 'pvalb'),
    ('5', 'sst'),
    ('5', 'vip'),
    # ('6', 'unknown'),
    ('6', 'ntsr1'),
    ('6', 'pvalb'),
    ('6', 'sst'),
    ('6', 'vip'),
]
mouse_types = OrderedDict([(typ, ('L'+typ[0], typ[1])) for typ in mouse_types])
mouse_types[('2/3', 'unknown')] = ('L2/3', 'pyr')

mouse_ee_types = OrderedDict([
    (('2/3', 'unknown'), 'L23pyr'),
    ((None, 'rorb'), 'rorb'),
    ((None, 'sim1'), 'sim1'),
    ((None, 'tlx3'), 'tlx3'),
    ((None, 'ntsr1'), 'ntsr1'),
])

mouse_nolayer_types = OrderedDict([
    ((None, 'unknown'), 'L23pyr'),
    ((None, 'cux2'), 'cux2'),
    ((None, 'rorb'), 'rorb'),
    ((None, 'nr5a1'), 'nr5a1'),
    ((None, 'sim1'), 'sim1'),
    ((None, 'tlx3'), 'tlx3'),
    ((None, 'rbp4'), 'rbp4'),
    ((None, 'slc17a8'), 'slc17a8'),
    ((None, 'ntsr1'), 'ntsr1'),
    ((None, 'ctgf'), 'ctgf'),
    ((None, 'pvalb'), 'pvalb'),
    ((None, 'sst'), 'sst'),
    ((None, 'vip'), 'vip'),
])

human_types = OrderedDict([
    (('1', 'unknown'), 'L1'),
    (('2', 'unknown'), 'L2'),
    (('3', 'unknown'), 'L3'), 
    (('4', 'unknown'), 'L4'), 
    (('5', 'unknown'), 'L5'), 
    (('6', 'unknown'), 'L6'),
])
human_types = OrderedDict([(typ, "L%s %s" % typ) for typ in human_types])

if args.organism == 'mouse':
    m0 = expts.matrix(mouse_types, mouse_types, mode='progress', title='Mouse Coarse Progress Matrix')
    m1 = expts.matrix(mouse_types, mouse_types, title='Mouse Coarse Cp Matrix')
    m2 = expts.matrix(mouse_ee_types, mouse_ee_types, title='Mouse L5 Recurrent Excitatory Matrix')
    m3 = expts.matrix(mouse_nolayer_types, mouse_nolayer_types, title='Mouse cre-type Matrix')

    # print total % complete
    matrix = expts.connectivity_matrix(mouse_types, mouse_types)
    n_started = (matrix['probed'] > 0).sum()
    probed = matrix['probed'] / 80.
    connected = matrix['connected'] / 10.
    progress = np.clip(np.where(probed > connected, probed, connected), 0, 1)
    print("Total progress %0.1f%%  (%d started)" % (100 * progress.sum()/90., n_started))

    #mouse distance plots
    plots = expts.distance_plot(pre_types=[('2/3', None)], post_types=[('2/3', None)], color=(255, 153, 51))
    expts.distance_plot('rorb', 'rorb', plots=plots, color=(102, 255, 102))
    expts.distance_plot('sim1', 'sim1', plots=plots, color=(102, 255, 255))
    expts.distance_plot('ntsr1', 'ntsr1', plots=plots, color=(153, 51, 255))

elif args.organism == 'human':
    m1 = expts.matrix(human_types, human_types, title='Human Cp Matrix')
    # m1 = expts.matrix(human_types, human_types, mode='progress')

    #human distance plots
    plots = expts.distance_plot(pre_types=[('2', None)], post_types=[('2', None)], color=(255, 153, 153))
    expts.distance_plot(pre_types=[('3', None)], post_types=[('3', None)], plots=plots, color=(255, 255, 102))
    # expts.distance_plot(pre_types=[('4', None)], post_types=[('4', None)], plots=plots, color=(102, 255, 102))
    expts.distance_plot(pre_types=[('5', None)], post_types=[('5', None)], plots=plots, color=(102, 255, 255))
    # expts.distance_plot(pre_types=[('6', None)], post_types=[('6', None)], plots=plots, color=(153, 51, 255))

