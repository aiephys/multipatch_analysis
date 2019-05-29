from __future__ import division, print_function

import os, sys, user, argparse
from multipatch_analysis import database as db
from multipatch_analysis.config import synphys_db
from multipatch_analysis import config
from multipatch_analysis.pipeline import all_modules

parser = argparse.ArgumentParser()
parser.add_argument('--reset-db', action='store_true', default=False, help="Drop all tables in the database.", dest='reset_db')
parser.add_argument('--vacuum', action='store_true', default=False, help="Ask the database to clean/optimize itself.")
parser.add_argument('--bake', type=str, default=None, help="Bake current database into an sqlite file.")
parser.add_argument('--skip-tables', type=str, default="", help="Comma-separated list of tables to skip while baking.", dest="skip_tables")
parser.add_argument('--overwrite', action='store_true', default=False, help="Overwrite existing sqlite file.")
parser.add_argument('--update', action='store_true', default=False, help="Update existing sqlite file.")
parser.add_argument('--dbg', action='store_true', default=False, help="Start debugging console.")

args = parser.parse_args(sys.argv[1:])

if args.dbg:
    import pyqtgraph as pg
    pg.dbg()

if args.reset_db:
    ans = raw_input('Reset database "%s"? ' % db.db_name)
    if ans == 'y':
        print("  Ok, here we go..")
        db.reset_db()
        print("    ..done.")
    else:
        print("  Oh very well. Some other time, perhaps.")

if args.vacuum:
    # cleans up DB and analyzes column statistics to improve query performance
    print("Mopping up %s.." % db.db_name)
    db.vacuum()
    print("   ..done.")


if args.bake is not None:
    if os.path.exists(args.bake):
        if args.overwrite:
            os.remove(config.synphys_db_sqlite)
        elif not args.update:
            print("sqlite database file %s already exists" % args.bake)
            sys.exit(0)
        
    db.bake_sqlite(args.bake, skip_tables=args.skip_tables.split(','))
