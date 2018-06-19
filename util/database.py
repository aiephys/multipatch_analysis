import sys, user
import multipatch_analysis.database.database as db
from multipatch_analysis.config import synphys_db

if '--reset-db' in sys.argv:
    ans = raw_input('Reset database "%s"? ' % db.db_name)
    if ans == 'y':
        print("  Ok, here we go..")
        db.reset_db()
        print("    ..done.")
    else:
        print("  Oh very well. Some other time, perhaps.")

if '--vacuum' in sys.argv:
    # cleans up DB and analyzes column statistics to improve query performance
    print("Mopping up %s.." % db.db_name)
    db.vacuum()
    print("   ..done.")
