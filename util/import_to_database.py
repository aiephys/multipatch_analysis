from __future__ import print_function

import os, sys, time, glob, argparse
import multiprocessing

import pyqtgraph as pg
pg.dbg()

from multipatch_analysis.experiment import Experiment
from multipatch_analysis.database.submission import SliceSubmission, ExperimentDBSubmission
from multipatch_analysis.database import database
from multipatch_analysis import config, synphys_cache, experiment_list, constants


def submit_expt(expt_id, raise_exc=False):
    # print(os.getpid(), expt_id, "start")
    global all_expts
    try:
        site_path = all_expts[expt_id]
        expt = Experiment(site_path=site_path)
        print("submit experiment: %0.3f" % expt_id, expt)
        
        slice_dir = expt.slice_dir
        sub = SliceSubmission(slice_dir)
        if not sub.submitted():
            sub.submit()
        
        start = time.time()
        
        sub = ExperimentDBSubmission(expt)
        if sub.submitted():
            print("   expt %s already in DB" % expt)
        else:
            sub.submit()
            print("   expt %s done (%s)" % (expt_id, expt))

        print("    %g sec" % (time.time()-start))
        return (expt_id, None)
    except Exception as exc:
        print(">>>> %d Error importing experiment %0.3f" % (os.getpid(), expt_id))
        sys.excepthook(*sys.exc_info())
        print("<<<< %0.3f" % expt_id)
        if raise_exc:
            raise
        return (expt_id, str(exc))
    # print(os.getpid(), expt_id, "return")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Import new experiments into the database.")
    parser.add_argument('--limit', type=int, default=None)
    parser.add_argument('--local', action='store_true', default=False)
    parser.add_argument('--workers', type=int, default=6)
    parser.add_argument('--uid', type=str, default=None)
    parser.add_argument('--raise-exc', action='store_true', default=False, dest='raise_exc', help='Do not ignore exceptions')
    
    args = parser.parse_args(sys.argv[1:])
    
    cache = synphys_cache.get_cache()
    all_expts = cache.list_experiments()

    if args.uid is not None:
        selected_expts = [float(uid) for uid in args.uid.split(',')]
    elif args.limit is not None:
        # Just a dirty trick to give a wider variety of experiments when we are testing
        # on a small subset
        import random
        random.seed(0)
        selected_expts = list(all_expts.keys())
        random.shuffle(selected_expts)
        selected_expts = selected_expts[:args.limit]
    else:
        selected_expts = list(all_expts.keys())

    print("Found %d cached experiments, will import %d." % 
          (len(all_expts), len(selected_expts)))
    print(selected_expts)
    
    if args.local is True:
        errors = []
        for i, expt in enumerate(selected_expts):
            errors.append(submit_expt(expt, raise_exc=args.raise_exc))
    else:
        ids = [expt for expt in selected_expts]

        # Dispose DB engine before forking, otherwise child processes will
        # inherit and muck with the same connections. See:
        # http://docs.sqlalchemy.org/en/rel_1_0/faq/connections.html#how-do-i-use-engines-connections-sessions-with-python-multiprocessing-or-os-fork
        database.engine.dispose()
        
        pool = multiprocessing.Pool(processes=args.workers, maxtasksperchild=1)
        errors = pool.map(submit_expt, ids, chunksize=1)  # note: maxtasksperchild is broken unless we also force chunksize

    errors = [e for e in errors if e[1] is not None]
    print("======= DB import complete with %d/%d errors =========" % (len(errors), len(selected_expts)))
    for expt_id, err in errors:
        print("%0.3f\t%s" % (expt_id, err))
    print("===================================================")
