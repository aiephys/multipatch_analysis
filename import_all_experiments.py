from __future__ import print_function

import os, sys, time, glob
from acq4.util.DataManager import getFileHandle
import pyqtgraph.multiprocess as mp

from experiment_list import ExperimentList
from submission import SliceSubmission, ExperimentDBSubmission
import database
import config
import synphys_cache


if __name__ == '__main__':
    import pyqtgraph as pg
    pg.dbg()
    
    print("Reading experiment list..")
    cache = synphys_cache.SynPhysCache()
    
    all_expts = cache.list_nwbs()
    
    for expt in all_expts:
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
        
    