from __future__ import print_function
import argparse, sys
import pyqtgraph as pg 
from multipatch_analysis.morphology import morphology_tables, init_tables, update_morphology
import multipatch_analysis.database as db


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Import morphological features into DB.")
    parser.add_argument('--rebuild', action='store_true', default=False, help="Remove and rebuild tables for this analysis")
    parser.add_argument('--workers', type=int, default=None, help="Set the number of concurrent processes during update")
    parser.add_argument('--local', action='store_true', default=False, help="Disable concurrent processing to make debugging easier")
    parser.add_argument('--raise-exc', action='store_true', default=False, help="Disable catching exceptions encountered during processing", dest='raise_exc')
    parser.add_argument('--limit', type=int, default=None, help="Limit the number of experiments to process")
    parser.add_argument('--expts', type=lambda s: [float(x) for x in s.split(',')], default=None, help="Select specific experiment IDs to analyze", )
    
    args = parser.parse_args(sys.argv[1:])
    if args.rebuild:
        args.rebuild = raw_input("Rebuild %s morphology table? " % db.db_name) == 'y'

    if args.local:
        pg.dbg()

    if args.rebuild:
        morphology_tables.drop_tables()
        init_tables()

    update_morphology(limit=args.limit, expts=args.expts, parallel=not args.local, workers=args.workers, raise_exceptions=args.raise_exc)
