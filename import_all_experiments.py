from __future__ import print_function

import os, sys, time, glob, argparse
from acq4.util.DataManager import getFileHandle
import pyqtgraph.multiprocess as mp
import multiprocessing

from experiment_list import ExperimentList
from submission import SliceSubmission, ExperimentDBSubmission
import database
import config
import synphys_cache


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--limit', type=int, default=None)
    parser.add_argument('--local', action='store_true', default=False)
    args, extra = parser.parse_known_args(sys.argv[1:])
    
    import pyqtgraph as pg
    pg.dbg()
    
    print("Reading experiment list..")
    cache = synphys_cache.SynPhysCache()
    
    all_nwbs = cache.list_nwbs()
    
    # Just a dirty trick to give a wider variety of experiments when we are testing
    # on a small subset
    import random
    random.seed(0)
    random.shuffle(all_nwbs)
    
    # clean out experiments with multiple NWB files for now
    all_sites = [os.path.dirname(nwb) for nwb in all_nwbs]
    n_sites = len(set(all_sites))
    all_sites = [expt for expt in all_sites if all_sites.count(expt) == 1]
    all_expts = [nwb for nwb in all_nwbs if os.path.dirname(nwb) in all_sites]
    
    if args.limit > 0:
        selected_expts = all_expts[:args.limit]
    else:
        selected_expts = all_expts
        
    print("Found %d nwb files in %d sites, will import %d experiments." % 
          (len(all_nwbs), n_sites, len(selected_expts)))
    
    for i, expt in enumerate(selected_expts):
        
        nwb_file = getFileHandle(expt)
        
        nwb_cache_file = getFileHandle(cache.get_cache(nwb_file.name()))
        site_dir = nwb_file.parent()
        slice_dir = site_dir.parent()

        try:
            start = time.time()
            
            
            print("submit slice:", slice_dir.name())
            sub = SliceSubmission(slice_dir)
            if sub.submitted():
                print("   already in DB")
            else:
                sub.submit()
            
            print("    %g sec" % (time.time()-start))
            start = time.time()
            
            print("submit site:", site_dir.name())
            print("    " + expt)
            sub = ExperimentDBSubmission(site_dir, nwb_cache_file)
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
                        rdm = proc._import('acq4.util.DataManager')
                        nwb_file = rdm.getFileHandle(nwb_file.name())
                        nwb_cache_file = rdm.getFileHandle(nwb_cache_file.name())
                        site_dir = nwb_file.parent()
                        slice_dir = site_dir.parent()
                        
                        rsub = proc._import('submission')
                        sub = rsub.ExperimentDBSubmission(site_dir, nwb_cache_file)
                        fn = sub.submit(_timeout=None)
                        
                    finally:
                        proc.join()

            print("    %g sec" % (time.time()-start))
        except Exception:
            print("=============================")
            print("failed:\n  %s" % (expt))
            #raise
            sys.excepthook(*sys.exc_info())
        
    