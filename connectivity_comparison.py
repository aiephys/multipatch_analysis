# *-* coding: utf-8 *-*

"""
Script comparing across multipatch experimental conditions.

"""

from __future__ import print_function, division

import argparse
import sys
import pyqtgraph as pg

from experiment_list import ExperimentList

parser = argparse.ArgumentParser()
parser.add_argument('--cre-type', type=str, help='Enter as pretype-posttype. If comparing connection types separate'
                    '' 'with ",". Ex pvalb-pavalb,pvalb-sst')
parser.add_argument('--calcium', action='store_true', default=False, dest='calcium',
                    help='cre-type must also be specified')
parser.add_argument('--age', type=str, help='Enter age ranges separated by ",". Ex 40-50,60-70.'
                    '' 'cre-type must also be specified')

args = parser.parse_args(sys.argv[1:])

cache_file = 'expts_cache.pkl'
all_expts = ExperimentList(cache=cache_file)
app = pg.mkQApp()

if args.cre_type is not None and len(args.cre_type.split(',')) == 1:
    cre_type = args.cre_type.split('-')
    if args.calcium is True:
        expts = all_expts.select(calcium='high')
        plots = expts.distance_plot(cre_type[0], cre_type[1], color=(0, 10), name=("%s->%s, calcium = 2.0mM " %(cre_type[0], cre_type[1])))
        expts = all_expts.select(calcium='low')
        expts.distance_plot(cre_type[0], cre_type[1], plots=plots, color=(5, 10), name=("%s->%s, calcium = 1.3mM " %(cre_type[0], cre_type[1])))
    elif args.age is not None:
        ages = args.age.split(',')
        if len(ages) < 2:
            print("Please specify more than one age range")
        expts = all_expts.select(age=ages[0])
        plots = expts.distance_plot(cre_type[0], cre_type[1], color=(0, 10), name=("%s->%s, age = P%s " %(cre_type[0], cre_type[1], ages[0])))
        expts = all_expts.select(age=ages[1])
        expts.distance_plot(cre_type[0], cre_type[1], plots=plots, color=(5, 10), name=("%s->%s, age = P%s " %(cre_type[0], cre_type[1], ages[1])))
elif args.cre_type is None and (args.calcium is not None or args.age is not None):
    print('Error: in order to compare across conditions a single cre-type connection must be specified')
else:
    cre_types = args.cre_type.split(',')
    plots = None
    for i, type in enumerate(cre_types):
        cre_type = type.split('-')
        expts = all_expts.select(cre_type=cre_type)
        plots = expts.distance_plot(cre_type[0], cre_type[1], plots=plots, color=(i, len(cre_types)*1.3))