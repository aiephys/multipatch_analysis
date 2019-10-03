import sys
import pyqtgraph as pg
from aisynphys.ui.dashboard import Dashboard

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--no-thread', action='store_true', default=False, dest='no_thread',
                    help='Do all polling in main thread (to make debugging easier).')
    parser.add_argument('--limit', type=int, dest='limit', default=0, help="Limit the number of experiments to poll (to make testing easier).")
    args = parser.parse_args(sys.argv[1:])

    app = pg.mkQApp()
    # console = pg.dbg()
    db = Dashboard(limit=args.limit, no_thread=args.no_thread)
    db.show()

    if sys.flags.interactive == 0:
        app.exec_()

