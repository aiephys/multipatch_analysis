from __future__ import print_function

import os, sys, time, glob, argparse
import multiprocessing

import pyqtgraph as pg
pg.dbg()

from multipatch_analysis.database.submission import SliceSubmission, ExperimentDBSubmission
from multipatch_analysis.database import database
from multipatch_analysis import config, synphys_cache, experiment_list


all_expts = experiment_list.cached_experiments()


def submit_expt(expt_id):
    try:
        expt = all_expts[expt_id]
        start = time.time()
        
        slice_dir = expt.slice_dir
        print("submit slice:", slice_dir)
        sub = SliceSubmission(slice_dir)
        if sub.submitted():
            print("   already in DB")
        else:
            sub.submit()
        
        print("    %g sec" % (time.time()-start))
        start = time.time()
        
        print("submit experiment:")
        print("    ", expt)
        sub = ExperimentDBSubmission(expt)
        if sub.submitted():
            print("   already in DB")
        else:
            sub.submit()

        print("    %g sec" % (time.time()-start))
    except Exception:
        print(">>>> %d Error importing experiment %s" % (os.getpid(), expt.uid))
        sys.excepthook(*sys.exc_info())
        print("<<<< %s" % expt.uid)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--limit', type=int, default=None)
    parser.add_argument('--local', action='store_true', default=False)
    parser.add_argument('--workers', type=int, default=6)
    parser.add_argument('--uid', type=str, default=None)
    args, extra = parser.parse_known_args(sys.argv[1:])
    
    if args.uid is not None:
        selected_expts = [all_expts[uid] for uid in args.uid.split(',')]
    else:
        # Just a dirty trick to give a wider variety of experiments when we are testing
        # on a small subset
        import random
        random.seed(0)
        shuffled = list(all_expts)
        random.shuffle(shuffled)
        
        if args.limit > 0:
            selected_expts = shuffled[:args.limit]
        else:
            selected_expts = shuffled
        
    print("Found %d cached experiments, will import %d." % 
          (len(all_expts), len(selected_expts)))
    print([ex.uid for ex in selected_expts])
    
    if args.local is True:
        for i, expt in enumerate(selected_expts):
            submit_expt(expt.uid)
    else:
        ids = [expt.uid for expt in selected_expts]

        # Dispose DB engine before forking, otherwise child processes will
        # inherit and muck with the same connections. See:
        # http://docs.sqlalchemy.org/en/rel_1_0/faq/connections.html#how-do-i-use-engines-connections-sessions-with-python-multiprocessing-or-os-fork
        database.engine.dispose()
        
        pool = multiprocessing.Pool(processes=args.workers, maxtasksperchild=1)
        pool.map(submit_expt, ids)

