from __future__ import print_function
import argparse, sys, os, logging
import six
from aisynphys.pipeline import all_pipelines
from aisynphys.database import default_db as db
from aisynphys import config


if __name__ == '__main__':
    logging.basicConfig(format="%(message)s")
    logging.getLogger('aisynphys').setLevel(logging.INFO)
    all_pipelines = all_pipelines()
    
    parser = argparse.ArgumentParser(description="Process analysis pipeline jobs")
    parser.add_argument('pipeline', type=str, help="The name of the pipeline to run: %s" % ', '.join(list(all_pipelines.keys())))
    parser.add_argument('modules', type=str, nargs='*', help="The name of the analysis module(s) to run")
    parser.add_argument('--update', action='store_true', default=False, help="Process any jobs that are ready to be updated")
    parser.add_argument('--retry', action='store_true', default=False, help="During update, retry processing jobs that previously failed (implies --update)")
    parser.add_argument('--force-update', action='store_true', default=False, help="During update, reprocess all available jobs regardless of status (allowed only with --limit or --uids)")
    parser.add_argument('--report', action='store_true', default=False, help="Print a report of pipeline status and errors", )
    parser.add_argument('--rebuild', action='store_true', default=False, help="Remove and rebuild tables for selected modules")
    parser.add_argument('--workers', type=int, default=None, help="Set the number of concurrent processes during update")
    parser.add_argument('--local', action='store_true', default=False, help="Disable concurrent processing to make debugging easier")
    parser.add_argument('--debug', action='store_true', default=False, help="Enable debugging features: disable parallel processing, raise exception on first error, open debugging gui")
    parser.add_argument('--limit', type=int, default=None, help="Limit the number of experiments to process")
    parser.add_argument('--uids', type=lambda s: s.split(','), default=None, help="Select specific IDs to analyze (or drop)", )
    parser.add_argument('--drop', action='store_true', default=False, help="Drop selected analysis results (do not run updates)", )
    parser.add_argument('--vacuum', action='store_true', default=False, help="Run VACUUM ANALYZE on the database to optimize its query planner", )
    parser.add_argument('--bake', action='store_true', default=False, help="Bake an sqlite file after the pipeline update completes", )
    parser.add_argument('--info', action='store_true', default=False, help="Display information about the selected pipeline", )
    
    
    args = parser.parse_args(sys.argv[1:])

    if args.debug:
        args.local = True
        import pyqtgraph as pg 
        logging.getLogger('aisynphys').setLevel(logging.DEBUG)
        pg.dbg()

    try:
        pipeline = all_pipelines[args.pipeline](database=db, config=config)
    except KeyError:
        raise Exception("Could not find pipeline named %s. Options are: %s"%(args.pipeline, str(all_pipelines.keys())))

    all_modules = pipeline.sorted_modules()
    
    if 'all' in args.modules:
        modules = list(all_modules.values())
    else:
        modules = []
        for mod in args.modules:
            try:
                if mod.startswith(':'):
                    i = list(all_modules.keys()).index(mod[1:])
                    modules.extend(list(all_modules.values())[:i+1])
                elif mod.endswith(':'):
                    i = list(all_modules.keys()).index(mod[:-1])
                    modules.extend(list(all_modules.values())[i:])
                else:
                    modules.append(all_modules[mod])
            except (KeyError, ValueError):
                print('Unknown analysis module "%s"; options are: %s' % (mod, list(all_modules.keys())))
                sys.exit(-1)

    # sort topologically
    modules = [m for m in list(all_modules.values()) if m in modules]

    if args.report:
        print("----------------------------------------------")
        print("Pipeline: %s   DB: %s" % (args.pipeline, str(db)))
        print("----------------------------------------------")
        print(pipeline.report(modules, job_ids=args.uids))
    
    if args.rebuild:
        mod_names = ', '.join([module.name for module in modules])
        if six.moves.input("Rebuild modules in %s: %s? (y/n) " % (str(db), mod_names)) != 'y':
            print("  Nuts.")
            sys.exit(-1)

    if args.bake and os.path.exists(config.synphys_db_sqlite):
        msg = "sqlite database file %s already exists; ok to overwrite? (y/n) " % config.synphys_db_sqlite
        ans = six.moves.input(msg)
        if ans == 'y':
            print("  Ok, you asked for it..")
            os.remove(config.synphys_db_sqlite)
        else:
            print("  Phooey.")
            args.bake = False

    if args.rebuild or args.drop:
        # this call takes care of logging also
        pipeline.drop(modules=modules, job_ids=args.uids)
        print("  done.")
 
    if args.update or args.rebuild or args.retry or args.force_update:
        if args.force_update and (args.limit is None) and (args.uids is None):
            print("Force-update permitted only with --uids or --limit options. Try rebuilding instead?")
            sys.exit(-1)
        report = []
        for module in modules:
            print("=============================================")
            result = module.update(job_ids=args.uids, retry_errors=args.retry, limit=args.limit, 
                                   parallel=not args.local, workers=args.workers, debug=args.debug,
                                   force=args.force_update)
            report.append((module, result))
            
        if args.vacuum:
            print("Starting vacuum..")
            db.vacuum()
            print("   ..finished vacuum.")

        print("\n================== Error Report ===========================")
        for module, result in report:
            print("------ %s : %d errors -------" % (module.name, result['n_errors']))
            for job, err in result['errors'].items():
                print("    %s : %s" % (job, err))
        
            
        print("\n================== Update Report ===========================")
        for module, result in report:
            print("{name:20s}  dropped: {n_dropped:6d}  updated: {n_updated:6d} ({n_retry:6d} retry)  errors: {n_errors:6d}".format(name=module.name, **result))

    if args.bake:
        print("\n================== Bake Sqlite ===========================")
        db.bake_sqlite(config.synphys_db_sqlite)

    if args.info:
        print("Pipeline:", args.pipeline)
        print("  {:20s}  :  ok  : err  :  dependencies".format('module'))
        print("  -----------------------------------------------------------")
        for module in all_modules.values():
            status = module.job_status()
            success = len([j for j in status.values() if j[0] is True])
            err = len(status) - success
            print("  {:20s}  : {:<5d}: {:<5d}:  {}".format(module.name, success, err, ', '.join([m.name for m in module.dependencies])))
