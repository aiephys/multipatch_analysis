from __future__ import print_function
import argparse, sys
import pyqtgraph as pg 
from multipatch_analysis.connection_strength import connection_strength_tables, init_tables, update_connectivity
import multipatch_analysis.database as db


if __name__ == '__main__':
    import user

    parser = argparse.ArgumentParser(description="Analyze connectivity and other properties of pairs, "
                                               "store to connection_strength table.")
    parser.add_argument('--update', action='store_true', default=False, help="Update tables with analysis from new experiments")
    parser.add_argument('--rebuild', action='store_true', default=False, help="Remove and rebuild tables for this analysis")
    parser.add_argument('--workers', type=int, default=6, help="Set the number of concurrent processes during update")
    parser.add_argument('--local', action='store_true', default=False, help="Disable concurrent processing to make debugging easier")
    parser.add_argument('--limit', type=int, default=0, help="Limit the number of experiments to process")
    
    args = parser.parse_args(sys.argv[1:])
    if args.rebuild:
        args.rebuild = raw_input("Rebuild %s connectivity table? " % db.db_name) == 'y'

    pg.dbg()

    if args.rebuild:
        args.update = True
        connection_strength_tables.drop_tables()

    init_tables()

    if args.update:
        update_connectivity()
