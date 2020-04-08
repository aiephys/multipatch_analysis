import re, argparse, sys
from datetime import datetime
from aisynphys.database import default_db as db

if __name__ == '__main__':

    args = sys.argv[1:]

    start = datetime.strptime(args[0], "%y%m%d")
    stop = datetime.strptime(args[1], "%y%m%d")

    tubes = db.query(db.PatchSeq.tube_id).all()
    n_tubes = 0

    for tube in tubes:
        name_check = re.match(r'P(M|T)S4_(?P<date>\d{6})_(?P<tube_id>\d{3})_A01', tube[0])
        if name_check is not None:
            date=name_check.group('date')
            date=datetime.strptime(date, "%y%m%d")
            if date >= start and date <= stop:
                n_tubes+=1

    print("Number of tubes from %s - %s: %d" % (start.date(), stop.date(), n_tubes))