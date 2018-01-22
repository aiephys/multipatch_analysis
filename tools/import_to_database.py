from __future__ import print_function

import os, sys, time, glob, argparse
import pyqtgraph.multiprocess as mp
import multiprocessing

from multipatch_analysis.database.submission import SliceSubmission, ExperimentDBSubmission
from multipatch_analysis.database import database
from multipatch_analysis import config, synphys_cache, experiment_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--limit', type=int, default=None)
    parser.add_argument('--local', action='store_true', default=False)
    args, extra = parser.parse_known_args(sys.argv[1:])
    
    import pyqtgraph as pg
    pg.dbg()
    
    all_expts = list(experiment_list.cached_experiments())
    
    # Just a dirty trick to give a wider variety of experiments when we are testing
    # on a small subset
    import random
    random.seed(0)
    random.shuffle(all_expts)
    
    if args.limit > 0:
        selected_expts = all_expts[:args.limit]
    else:
        selected_expts = all_expts
        
    print("Found %d cached experiments, will import %d." % 
          (len(all_expts), len(selected_expts)))
    
    for i, expt in enumerate(selected_expts):
        try:
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

                if '--local' in sys.argv:
                    sub.submit()
                    
                else:
                    # Do experiment submission in a subprocess to control memory
                    # usage.
                    proc = mp.Process(pyqtapis={'QString': 2, 'QVariant': 2})
                    try:
                        rxl = proc._import('multipatch_analysis.experiment_list')
                        all_expts = rxl.cached_experiments()
                        expt = all_expts[expt.uid]

                        rsub = proc._import('multipatch_analysis.database.submission')
                        sub = rsub.ExperimentDBSubmission(expt)
                        fn = sub.submit(_timeout=None)
                        
                    finally:
                        proc.join()

            print("    %g sec" % (time.time()-start))
        except Exception:
            print("=============================")
            print("failed:\n  %s" % (expt))
            #raise
            sys.excepthook(*sys.exc_info())
        
    