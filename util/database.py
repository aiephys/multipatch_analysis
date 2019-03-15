import sys, user, argparse
from multipatch_analysis import database as db
from multipatch_analysis.config import synphys_db
from multipatch_analysis import config


parser = argparse.ArgumentParser()
parser.add_argument('--reset-db', action='store_true', default=False, help="Drop all tables in the database.", dest='reset_db')
parser.add_argument('--vacuum', action='store_true', default=False, help="Ask the database to clean/optimize itself.")
parser.add_argument('--bake', action='store_true', default=False, help="Bake current database into an sqlite file.")

args = parser.parse_args(sys.argv[1:])

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


if args.bake:
    sqlite_addr = "sqlite:///%s" % config.synphys_db_sqlite
    sqlite_engine = db.database.create_engine(sqlite_addr)
    db.database.create_tables(engine=sqlite_engine)
    