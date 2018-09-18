from __future__ import print_function
import argparse, sys
# sys.path.append('/home/corinnet/workspace/aiephys/multipatch_analysis/analyses')
# import strength_analysis
import pyqtgraph as pg 
from multipatch_analysis.first_pulse_fits import first_pulse_fit_tables, init_tables, update_fit
import multipatch_analysis.database as db


if __name__ == '__main__':
    import user

    parser = argparse.ArgumentParser(description="Fit individual first pulse psps "
                                               "store to single_pulse_fit tables.")
    parser.add_argument('--rebuild', action='store_true', default=False, help="Remove and rebuild tables for this analysis")
    parser.add_argument('--workers', type=int, default=6, help="Set the number of concurrent processes during update")
    parser.add_argument('--local', action='store_true', default=True, help="Disable concurrent processing to make debugging easier")
    parser.add_argument('--raise-exc', action='store_true', default=False, help="Disable catching exceptions encountered during processing", dest='raise_exc')
    parser.add_argument('--limit', type=int, default=None, help="Limit the number of experiments to process")
    
    args = parser.parse_args(sys.argv[1:])
    if args.rebuild:
        args.rebuild = raw_input("Rebuild '%s' first pulse fit tables? " % db.db_name) == 'y'

    pg.dbg()

    if args.rebuild:
        first_pulse_fit_tables.drop_tables()
        init_tables()

    update_fit(limit=100, expts=None, parallel=False, workers=6, raise_exceptions=False, session=None)
