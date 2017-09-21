from __future__ import print_function

import sys, time
from acq4.util.DataManager import getDirHandle, getFileHandle
from pyqtgraph.multiprocess.parallelizer import Parallelize

from submission import SliceSubmission, ExperimentDBSubmission
from experiment_list import ExperimentList



if __name__ == '__main__':
    import pyqtgraph as pg
    pg.dbg()
    
    cache_file = 'expts_cache.pkl'
    all_expts = ExperimentList(cache=cache_file)

    n = 0
    for expt in all_expts:
        nwb_file = getFileHandle(expt.nwb_file)
        nwb_cache_file = getFileHandle(expt.nwb_cache_file)
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
            sub = ExperimentDBSubmission(site_dir, nwb_cache_file)
            if sub.submitted():
                print("   already in DB")
            else:
                sub.submit()
                n += 1

            print("    %g sec" % (time.time()-start))
        except Exception:
            print("=============================")
            print("failed:\n  %s\n  %s" % (slice_dir, site_dir))
            #raise
            sys.excepthook(*sys.exc_info())
            
        if n > 15:
            sys.exit("that's enough for now.")
        
    #tasks = [expt.uid for expt in all_expts]
    #with Parallelize(tasks, workers=4, results=results) as tasker:
        #for task in tasker:
            #result = processTask(task)
            #tasker.results.append(result)

    