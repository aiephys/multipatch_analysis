from multipatch_analysis.connection_strength import (
    connection_strength_tables, pulse_response_strength_tables, init_tables, rebuild_strength, rebuild_connectivity
)


if __name__ == '__main__':
    import user

    parser = argparse.ArgumentParser()
    parser.add_argument('--rebuild', action='store_true', default=False, "Rebuild connection strength and connectivity tables")
    parser.add_argument('--rebuild-connectivity', action='store_true', default=False, dest='rebuild_connectivity', "Rebuild connectivity table, leaving connection strength table unaffected")
    parser.add_argument('--local', action='store_true', default=False, "Disable concurrent processing to make debugging easier")
    parser.add_argument('--workers', type=int, default=6, help="Set the number of concurrent processes during rebuild")
    parser.add_argument('--seed', type=int, default=-1, help="Random seed used to shuffle classifier training data")
    parser.add_argument('--limit', type=int, default=0, help="Limit the number of experiments to process")
    
    args = parser.parse_args(sys.argv[1:])
    if args.rebuild:
        args.rebuild = raw_input("Rebuild %s strength tables? " % db.db_name) == 'y'
    if args.rebuild_connectivity:
        args.rebuild_connectivity = raw_input("Rebuild %s connectivity table? " % db.db_name) == 'y'

    pg.dbg()

    if args.rebuild:
        connection_strength_tables.drop_tables()
        pulse_response_strength_tables.drop_tables()
        init_tables()
        rebuild_strength(limit=limit, parallel=(not args.local), workers=args.workers)
        rebuild_connectivity()
    elif args.rebuild_connectivity:
        print("drop tables..")
        connection_strength_tables.drop_tables()
        print("create tables..")
        init_tables()
        print("rebuild..")
        rebuild_connectivity()
    else:
        init_tables()

