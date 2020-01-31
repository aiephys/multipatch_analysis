import re, argparse, sys
from datetime import datetime
from aisynphys.database import default_db as db

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--month', nargs=1, type=int)
    parser.add_argument('--year', nargs=1, type=int)

    args = parser.parse_args(sys.argv[1:])

    month = args.month
    year = args.year

    if year is None:
        year = 2020
    else:
        year = year[0]

    tubes = db.query(db.PatchSeq.tube_id).all()
    n_tubes = 0

    for tube in tubes:
        name_check = re.match(r'P(M|T)S4_(?P<date>\d{6})_(?P<tube_id>\d{3})_A01', tube[0])
        if name_check is not None:
            date=name_check.group('date')
            date=datetime.strptime(date, "%y%m%d")
            if month is not None:
                if date.month==month[0] and date.year==year:
                    n_tubes+=1
            else:
                if date.year == year:
                    n_tubes+=1

    print("Number of tubes for %s %s: %d" % (month, year, n_tubes))