import sys, multiprocessing
import numpy as np
from collections import OrderedDict
from pyqtgraph import toposort


class PipelineModule(object):
    """Pipeline modules represent analysis tasks that can be run independently of other parts of the analysis
    pipeline. 
    
    For any given experiment, a sequence of analysis stages must be processed in order. Each stage requires a specific
    set of inputs to be present, and produces/stores some output. Inputs and outputs can be raw data files, database tables,
    etc., although outputs are probably always written into a database. Each PipelineModule subclass represents a single stage
    in the sequence of analyses that occur across a pipeline.
    
    The work done by a single stage in the pipeline is divided up into jobs (units of work), where each stage may
    decide for itself what a suitable unit of work is. For many stages, the unit of work is the analysis done on a single experiment.
    Some stages may have finer granularity (multiple jobs per experiment), and some stages may have coarser granularity
    (multiple experiments per job), or even only have a single unit of work for the entire database, such as when aggregating 
    very high-level results.
    
    Note that PipelineModule classes generally need not implement any actual _analysis_; rather, they are simply responsible for
    data _management_ within the analysis pipeline. When an PipelineModule is asked to update an analysis result, it may call out to
    other packages to do the real work.
    
    From here, we should be able to:
    - Ask for which jobs this analysis has already been run (and when)
    - Ask for which jobs this analysis needs to be run (input dependencies are met and no result exists yet or is out of date)    
    - Run and store analysis for any specific experiment
    - Remove / re-run analysis for any specific experiment
    - Remove / re-run analysis for all jobs (for example, after a code change affecting all analysis results)
    - Initialize storage for analysis results (create DB tables, empty folders, etc.)
    
    """    
    
    name = None
    dependencies = []

    @staticmethod
    def all_modules():
        subclasses = PipelineModule.__subclasses__()
        deps = {c:c.dependencies for c in subclasses}
        return OrderedDict([(mod.name, mod) for mod in toposort(deps)])
    
    @classmethod
    def dependent_modules(cls):
        mods = cls.all_modules()
        return [mod for mod in mods if cls in mod.dependencies]
    
    @classmethod
    def update(cls, job_ids=None, limit=None, parallel=False, workers=None, raise_exceptions=False):
        """Update analysis results for this module.
        
        Parameters
        ----------
        job_ids : list | None
            List of job IDs to be updated, or None to update all jobs.
        parallel : bool
            If True, run jobs in parallel threads or subprocesses.
        workers : int or None
            Number of parallel workers to spawn. If None, then use one worker per CPU core.
        limit : int | None
            Maximum number of jobs to process (or None to disable this limit).
            If limit is enabled, then jobs are randomly shuffled before selecting the limited subset.
        raise_exceptions : bool
            If True, then exceptions are raised and will end any further processing.
            If False, then errors are logged and ignored.
            This is used mainly for debugging to allow traceback inspection.
        """
        print("Updating pipeline stage: %s" % cls.name)
        if job_ids is None:
            print("Searching for jobs to update..")
            run_job_ids, drop_job_ids = cls.updatable_jobs()                    
            if limit is not None:
                np.random.shuffle(job_ids)
                job_ids = job_ids[:limit]
        else:
            run_job_ids = job_ids
            drop_job_ids = job_ids
            
        print("Found %d jobs to update." % len(run_job_ids))

        if len(drop_job_ids) > 0:
            print("Dropping %d invalid results.." % len(drop_job_ids))
            cls.drop_jobs(drop_job_ids)

        run_jobs = [(expt_id, i, len(run_job_ids)) for i, expt_id in enumerate(run_job_ids)]

        if parallel:
            print("Processing all jobs (parallel)..")
            pool = multiprocessing.Pool(processes=workers)
            pool.map(cls._run_job, run_jobs)
        else:
            print("Processing all jobs (serial)..")
            for job in run_jobs:
                cls._run_job(job, raise_exceptions=raise_exceptions)

    @classmethod
    def _run_job(cls, job, raise_exceptions=False):
        """Entry point for running a single analysis job; may be invoked in a subprocess.
        """
        job_id, job_index, n_jobs = job
        print("Processing %d/%d  %s") % (job_index, n_jobs, job_id)
        try:
            cls.process_job(job_id)
        except Exception:
            if raise_exceptions:
                raise
            else:
                print("Error processing %d/%d  %s:") % (job_index, n_jobs, job_id)
                sys.excepthook(*sys.exc_info())
        else:
            print("Finished %d/%d  %s") % (job_index, n_jobs, job_id)
    
    @classmethod
    def process_job(cls, job_id):
        """Process analysis for one job.
        
        Parameters
        ----------
        job_id :
            The ID of the job to be processed.
        
        Must be reimplemented in subclasses.
        """
        raise NotImplementedError()

    @classmethod
    def initialize(cls):
        """Create space (folders, tables, etc.) for this analyzer to store its results.
        """
        raise NotImplementedError()
        
    @classmethod
    def drop_jobs(cls, job_ids):
        """Remove all results previously stored for a list of job IDs.
        """
        raise NotImplementedError()

    @classmethod
    def drop_all(cls):
        """Remove all results generated by this module.
        """
        raise NotImplementedError()

    @classmethod
    def finished_jobs(cls):
        """Return an ordered dict of job IDs that have been processed by this module and
        the dates when they were processed.

        Note that some results returned may be obsolete if dependencies have changed.
        """
        raise NotImplementedError()

    # @classmethod
    # def job_summary(cls):
    #     """Return information about jobs that have (not) been processed.
        
    #     Returns
    #     -------
    #     finished : list
    #         List of job IDs that have finished, valid results
    #     invalid : list
    #         List of job IDs that have finished, invalid results
    #     ready : list
    #         List of job IDs that have no result but are ready to be processed
    #     """
    #     my_finished = cls.finished_jobs()
    #     dep_finished = cls.ready_jobs()
        
    #     finished = []
    #     invalid = []
    #     ready = []
        
    #     all_jobs = sorted(list(set(list(my_finished.keys()) + list(dep_finished.keys()))))
        
    #     for job in all_jobs:
    #         if job in my_finished:
    #             if job not in dep_finished or dep_finished[job] > my_finished[job]:
    #                 invalid.append(job)
    #             else:
    #                 finished.append(job)
    #         else:
    #             ready.append(job)
                
    #     return finished, invalid, ready

    @classmethod
    def ready_jobs(cls):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        jobs = OrderedDict()
        for i,dep in enumerate(cls.dependencies):
            dep_finished = dep.finished_jobs()
            for k,v in dep_finished.items():
                jobs.setdefault(k, []).append(v)
        
        finished = OrderedDict()
        for job, dates in jobs.items():
            if len(dates) < len(cls.dependencies):
                # not all dependencies are finished yet
                continue
            finished[job] = max(dates)

        return finished

    @classmethod
    def updatable_jobs(cls):
        """Return lists of jobs that should be updated and/or should have their results dropped.
        """
        run_job_ids = []
        drop_job_ids = []
        ready = cls.ready_jobs()
        finished = cls.finished_jobs()
        for job in ready:
            if job in finished:
                if ready[job] > finished[job]:
                    # result is invalid
                    run_job_ids.append(job)
                    drop_job_ids
                else:
                    # result is valid
                    pass  
            else:
                # no current result
                run_job_ids.append(job)
        
        # look for orphaned results
        for job in finished:
            if job not in ready and job not in drop_job_ids:
                drop_job_ids.append(job)
        
        return run_job_ids, drop_job_ids