from __future__ import division, print_function

import os, sys, argparse
import six
from aisynphys.database import default_db as db
from aisynphys.config import synphys_db
from aisynphys import config

parser = argparse.ArgumentParser()
parser.add_argument('--reset-db', action='store_true', default=False, help="Drop all tables in the database.", dest='reset_db')
parser.add_argument('--vacuum', action='store_true', default=False, help="Ask the database to clean/optimize itself.")
parser.add_argument('--bake', type=str, default=None, help="Bake current database into an sqlite file.")
parser.add_argument('--clone', type=str, default=None, help="Clone current database into a new database with the given name.")
parser.add_argument('--tables', type=str, default=None, help="Comma-separated list of tables to include while baking.")
parser.add_argument('--skip-tables', type=str, default="", help="Comma-separated list of tables to skip while baking.", dest="skip_tables")
parser.add_argument('--skip-columns', type=str, default="", help="Comma-separated list of table.column names to skip while baking.", dest="skip_columns")
parser.add_argument('--overwrite', action='store_true', default=False, help="Overwrite existing sqlite file.")
parser.add_argument('--update', action='store_true', default=False, help="Update existing sqlite file.")
parser.add_argument('--drop', type=str, default=None, help="Drop database with the given name.")
parser.add_argument('--dbg', action='store_true', default=False, help="Start debugging console.")


args = parser.parse_args(sys.argv[1:])

if args.dbg:
    import pyqtgraph as pg
    pg.dbg()

if args.reset_db:
    ans = six.moves.input('Reset database "%s"? ' % db)
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


tables = None if args.tables is None else args.tables.split(',')
if args.bake is not None:
    if os.path.exists(args.bake):
        if args.overwrite:
            os.remove(args.bake)
        elif not args.update:
            print("sqlite database file %s already exists" % args.bake)
            sys.exit(0)
        
    skip_cols = {}
    for colname in args.skip_columns.split(','):
        table, col = colname.split('.')
        skip_cols.setdefault(table, []).append(col)
    db.bake_sqlite(args.bake, tables=tables, skip_tables=args.skip_tables.split(','), skip_columns=skip_cols)


if args.clone is not None:
    db.clone_database(args.clone, tables=tables, skip_tables=args.skip_tables.split(','))


if args.drop is not None:
    drop_db = db.get_database(args.drop)
    ans = six.moves.input("Seriously? I'm gonna dump the entire \"%s\" database? (y/n) " % drop_db)
    if ans == 'y':
        print("  Ok, don't look..")
        drop_db.drop_database()
        print("    ..done.")
    else:
        print("  Whew! Close one.")
    
