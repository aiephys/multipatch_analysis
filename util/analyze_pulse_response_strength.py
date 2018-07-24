from __future__ import print_function
import argparse, sys
import pyqtgraph as pg 
from multipatch_analysis.pulse_response_strength import pulse_response_strength_tables, init_tables, update_strength
import multipatch_analysis.database as db


if __name__ == '__main__':
    import user

    parser = argparse.ArgumentParser(description="Analyze strength and other properties of individual pulse responses, "
                                               "store to pulse_response_strength and baseline_response_strength tables.")
    parser.add_argument('--rebuild', action='store_true', default=False, help="Remove and rebuild tables for this analysis")
    parser.add_argument('--workers', type=int, default=6, help="Set the number of concurrent processes during update")
    parser.add_argument('--local', action='store_true', default=False, help="Disable concurrent processing to make debugging easier")
    parser.add_argument('--raise-exc', action='store_true', default=False, help="Disable catching exceptions encountered during processing", dest='raise_exc')
    parser.add_argument('--limit', type=int, default=0, help="Limit the number of experiments to process")
    parser.add_argument('--expts', type=lambda s: [float(x) for x in s.split(',')], default=None, help="Select specific experiment IDs to analyze", )
    
    args = parser.parse_args(sys.argv[1:])
    if args.rebuild:
        args.rebuild = raw_input("Rebuild %s pulse response strength tables? " % db.db_name) == 'y'

    if args.local:
        pg.dbg()

    if args.rebuild:
        pulse_response_strength_tables.drop_tables()
    
    init_tables()

    update_strength(limit=args.limit, expts=args.expts, parallel=(not args.local), workers=args.workers, raise_exceptions=args.raise_exc)
