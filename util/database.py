import os, sys, user, argparse
from multipatch_analysis import database as db
from multipatch_analysis.config import synphys_db
from multipatch_analysis import config
from multipatch_analysis.pipeline import all_modules

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
    if os.path.exists(config.synphys_db_sqlite):
        msg = "sqlite database file %s already exists; ok to overwrite? " % config.synphys_db_sqlite
        ans = raw_input(msg)
        if ans == 'y':
            print("  Ok, you asked for it..")
            os.remove(config.synphys_db_sqlite)
        else:
            print("  Phooey.")
            sys.exit(0)
        
    sqlite_addr = "sqlite:///%s" % config.synphys_db_sqlite
    sqlite_engine = db.database.create_engine(sqlite_addr)
    db.database.create_tables(engine=sqlite_engine)
    
    read_session = db.Session()
    write_session = db.database.sessionmaker(bind=sqlite_engine)()
    for mod in all_modules().values():
        table_group = mod.table_group
        for table in table_group.schemas:
            print("Querying %s.." % table)
            recs = read_session.query(table_group[table]).all()
            print("   pulled %d records, writing.." % len(recs))
            for rec in recs:
                write_session.merge(rec, load=False)
            write_session.commit()
            print("   done with %s." % table)
    
    print("All finished!")
            
            